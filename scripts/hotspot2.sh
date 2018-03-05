#! /bin/bash
set -x -e -o pipefail

usage() {
  cat >&2 <<__EOF__
Usage:  $0 [options] in.bam outdir

Options:
    -h                    Show this helpful help

  Mandatory options:
    -c CHROM_SIZES_FILE   (lowercase 'c')
                          File containing the lengths of all chromosomes
                          to include in the analysis, in BED or starch format.
                          All start coordinates (column 2) must be 0.
    -C CENTER_SITES_FILE  (uppercase 'C')
                          File of mappble sites (1bp each) where cleavages
                          observed at mappable sites within a radius of
                          NEIGHBORHOOD_SIZE bp will be tallied
                          and used to call hotspots of cleavage activity.
                          IMPORTANT: The user needs to create this file
                          using the script extractCenterSites.sh before running
                          $0, and the same NEIGHBORHOOD_SIZE
                          must be specified for both scripts.

  Optional options (note distinction between 'f' and 'F'):
    -M MAPPABLE_REG_FILE  (uppercase 'M')
                          The file of mappable regions that was used
                          to create the CENTER_SITES_FILE.
    -n NEIGHBORHOOD_SIZE  Local neighborhood size (100)
    -w WINDOW_SIZE        Background region size  (25000)
    -m MIN_HOTSPOT_WIDTH  Minimum hotspot width allowed (50)
    -P                    Write P-values to output file, in column 6 (not written by default)
    -f HOTSPOT_THRESHOLD  False-discovery rate to use for hotspot filtering (0.05)
    -F SITECALL_THRESHOLD False-discovery rate to use for site-call filtering (0.05)
    -S SMOOTHING_PARAM    (uppercase 'S')
                          Advanced option, to influence curve fitting (5).
                          Should be a small odd integer >=5; for noisy data, try 17.

    Neighborhood and window sizes are specified as the distance from the edge
    to the center - i.e, a 100bp neighborhood size is a 201bp window.

    The site-call False Discovery Rate threshold (-F) must be greater than or
    equal to the hotspot FDR threshold (-f).  After successful completion of this script,
    the user can, if desired, efficiently call hotspots at any threshold value
    HOTSPOT_THRESHOLD < x <= SITECALL_THRESHOLD without re-running $0,
    via the script hsmerge.sh. It is generally recommended to set SITECALL_THRESHOLD
    to the lowest value at which you might to investigate hotspots, e.g. 0.05 or 0.10.
    Higher values are generally only useful for debugging.

__EOF__
  exit 2
}

log() {
  echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t$*"
}

require_exes() {
  for x in "$@"; do
    if ! which "$x" &>/dev/null; then
      echo "Could not find $x!"
      exit -1
    fi
  done
}

CHROM_SIZES=""
CENTER_SITES=""
MAPPABLE_REGIONS=""
SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=100 # i.e., 201bp regions
BACKGROUND_WINDOW_SIZE=50001           # i.e., +/-25kb around each position
MIN_HOTSPOT_WIDTH=50
HOTSPOT_FDR_THRESHOLD="0.05"
CALL_THRESHOLD="$HOTSPOT_FDR_THRESHOLD"
SEED=$RANDOM
WRITE_PVALS=""
SMOOTHING_PARAM=""

# Note: Options in the string that are not immediately followed by ':'
# will not get a value read for them.  Examples are h and P.
while getopts 'hc:C:M:e:f:F:m:n:p:s:S:w:P' opt; do
  case "$opt" in
    h)
      usage
      ;;
    c)
      CHROM_SIZES=$OPTARG
      ;;
    C)
      CENTER_SITES=$OPTARG
      ;;
    M)
      MAPPABLE_REGIONS=$OPTARG
      ;;
    f)
      HOTSPOT_FDR_THRESHOLD=$OPTARG
      ;;
    F)
      CALL_THRESHOLD=$OPTARG
      ;;
    m)
      MIN_HOTSPOT_WIDTH=$OPTARG
      ;;
    n)
      SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=$OPTARG
      ;;
    P)
      WRITE_PVALS="--write_pvals"
      ;;
    s)
      SEED=$OPTARG
      ;;
    S)
      SMOOTHING_PARAM="-m $OPTARG"
      ;;
    w)
      BACKGROUND_WINDOW_SIZE=$((2 * OPTARG + 1))
      ;;

  esac
done
shift $((OPTIND - 1))

if [ "$CHROM_SIZES" == "" ]; then
  echo -e "Error:  Required argument -c CHROM_SIZES_FILE was not provided."
  usage
fi

if [ ! -s "$CHROM_SIZES" ]; then
  echo -e "Error:  CHROM_SIZES file \"$CHROM_SIZES\" was not found, or is empty."
  usage
fi

if [ "$CENTER_SITES" == "" ]; then
  echo -e "Error:  Required argument -C CENTER_SITES_FILE was not provided."
  usage
fi

if [ ! -s "$CENTER_SITES" ]; then
  echo -e "Error:  CENTER_SITES file \"$CENTER_SITES\" was not found, or is empty."
  usage
fi

if [ "$MAPPABLE_REGIONS" != "" ] && [ ! -s "$MAPPABLE_REGIONS" ]; then
  echo -e "Error:  MAPPABLE_REGIONS file \"$MAPPABLE_REGIONS\" was not found, or is empty."
  usage
fi

# Check to make sure Hotspot FDR <= site-call FDR
if awk '{exit $1>$2?0:1}' <<<"$HOTSPOT_FDR_THRESHOLD $CALL_THRESHOLD"; then
  echo "Hotspot FDR threshold (-f $HOTSPOT_FDR_THRESHOLD) cannot be greater than site-calling threshold (-F $CALL_THRESHOLD)" >&2
  exit 2
fi

if [[ -z "$1" || -z "$2" ]]; then
  usage
fi

BAM=$1
OUTDIR=$2

log "Checking system for required executables..."
require_exes modwt bedGraphToBigWig bedmap samtools hotspot2_part1 hotspot2_part2

WAVELETS_EXE=$(which modwt)
CUTCOUNT_EXE="$(dirname "$0")/cutcounts.bash"
DENSPK_EXE="$(dirname "$0")/density-peaks.bash"
MERGE_EXE="$(dirname "$0")/hsmerge.sh"
HOTSPOT_EXE1=hotspot2_part1
HOTSPOT_EXE2=hotspot2_part2

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)

mkdir -p "$OUTDIR"

base="$OUTDIR/$(basename "$BAM" .bam)"

HOTSPOT_OUTFILE="$base.hotspots.fdr$HOTSPOT_FDR_THRESHOLD.starch"
CUTCOUNTS="$base.cutcounts.starch"
FRAGMENTS_OUTFILE="$base.fragments.sorted.starch"
TOTALCUTS_OUTFILE="$base.cleavage.total"
OUTFILE="$base.allcalls.starch"
DENSITY_OUTFILE="$base.density.starch"
DENSITY_BW="$base.density.bw"
PEAKS_OUTFILE="$base.peaks.starch"
SPOT_SCORE_OUTFILE="$base.SPOT.txt"

clean=0
if [[ -z "$TMPDIR" ]]; then
  TMPDIR=$(mktemp -d)
  clean=1
fi

# temporary files
TEMP_CHROM_MAPPING_HOTSPOT2PART1=${TMPDIR}/temp_chrom_mapping_hotspot2part1.txt
TEMP_PVALS=${TMPDIR}/temp_pvals.txt
TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1=${TMPDIR}/temp_intermediateFile_hotspot2part1.txt

log "Generating cut counts..."
bash "$CUTCOUNT_EXE" "$BAM" "$CUTCOUNTS" "$FRAGMENTS_OUTFILE" "$TOTALCUTS_OUTFILE" "$CHROM_SIZES" "$MAPPABLE_REGIONS"

log "Tallying filtered cut counts in small windows and running part 1 of hotspot2..."
bedmap --faster --range "$SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE" --delim "\t" --prec 0 --echo --sum "$CENTER_SITES" "$CUTCOUNTS" \
  | "$AWK_EXE" 'BEGIN{OFS="\t"}{if("NAN"==$4){$4=0} print $1, $2, $3, "i", $4}' \
  | "$HOTSPOT_EXE1" --background_size="$BACKGROUND_WINDOW_SIZE" -c $TEMP_CHROM_MAPPING_HOTSPOT2PART1 -p $TEMP_PVALS $SMOOTHING_PARAM -o $TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1
if [ "$?" != "0" ]; then
    echo -e "An error occurred when calling bedmap on the \"center sites\" and filtered cut counts file, or while running part 1 of hotspot2."
    exit 2
fi

numEntries=`wc -l < $TEMP_PVALS` # used to aid memory allocation

log "Running part 2 of hotspot2..."
"$HOTSPOT_EXE2" --fdr_threshold="$CALL_THRESHOLD" $WRITE_PVALS \
		-i $TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1 -n $numEntries -c $TEMP_CHROM_MAPPING_HOTSPOT2PART1 -p $TEMP_PVALS \
    | starch - \
    >"$OUTFILE"

# We report the largest -log10(FDR) observed at any bp of a hotspot
# as the "score" of that hotspot, where FDR is the site-specific FDR estimate.
# FDR values of 0 will be encountered, and we don't want to do log(0).
# Nonzero FDR values as low as 1e-308 have been seen during testing.
# We choose to cap all FDR estimates at 1e-100, i.e. -log10(FDR) = 100.
# The constant c below converts from natural logarithm to log10.
# To combat issues awk sometimes has with numbers below 1e-300,
# we parse the FDR initially as a string, and assume anything below 1e-100
# passes the user's threshold.

log "Calling hotspots..."

"$MERGE_EXE" \
  -f "$HOTSPOT_FDR_THRESHOLD" \
  -m "$MIN_HOTSPOT_WIDTH" \
  "$OUTFILE" \
  "$HOTSPOT_OUTFILE"

log "Calculating SPOT score..."
num_cleaves=$(cat "$TOTALCUTS_OUTFILE")
cleaves_in_hotspots=$(bedops --ec -e 1 "$CUTCOUNTS" "$HOTSPOT_OUTFILE" | awk 'BEGIN{s=0} {s+=$5} END {print s}')
echo "scale=4; $cleaves_in_hotspots / $num_cleaves" \
  | bc \
    >"$SPOT_SCORE_OUTFILE"

#log "Creating peaks and density..."
bash "$DENSPK_EXE" "$TMPDIR" "$WAVELETS_EXE" "$CUTCOUNTS" "$HOTSPOT_OUTFILE" "$CHROM_SIZES" "$DENSITY_OUTFILE" "$PEAKS_OUTFILE"

log "Converting density to bigwig..."
TMPFRAGS="$(mktemp -t fragsXXXXX)"
unstarch "$DENSITY_OUTFILE" | cut -f1,2,3,5 >"$TMPFRAGS"
bedGraphToBigWig \
  "$TMPFRAGS" \
  <(cut -f1,3 "$CHROM_SIZES") \
  "$DENSITY_BW"

log "Done!"

rm -f "$TMPFRAGS"
if [[ $clean != 0 ]]; then
  rm -rf "$TMPDIR"
fi

exit 0
