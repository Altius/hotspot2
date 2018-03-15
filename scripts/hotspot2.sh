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
    -p PEAKS_DEFINITION   (lowercase 'p'; default = "default_peaks")
                          Specifies the definition of "peaks" (subsets of hotspots)
                          that will be used for the accompanying peaks output file.
                          If supplied, its value must be one of the following:
                          "default_peaks", "always_summit_centered", "varWidth_nn_ID",
                          with these exact capitalizations (or lack thereof).
                          "always_summit_centered" will return 150-bp fixed-width
                          peaks always centered on the summit of the corresponding
                          wavelet (or of the max density, when a hotspot lacks a wavelet summit).
                          Such peaks will sometimes extend beyond the edge of a hotspot.
                          "default_peaks" behaves similarly to "always_summit_centered",
                          except that 150-bp fixed-width peaks extending beyond the edge
                          of their parent hotspot are nudged to lie fully within the hotspot
                          whenever possible.
                          "varWidth_nn_ID" specifies that variable-width peaks,
                          down to a minimum width of nn bp, are to be reported.
                          (Use "varWidth_50_ID" for a minimum width of 50 bp, etc.)
                          "ID" is an identifier written into column 4 of the peaks output file.
                          These peaks correspond to full-width-at-half-maximum
                          elements for local maxima within the density curve;
                          a local minimum is substituted whenever half-maximum is not attained.
                          These variable-width peaks are guided by wavelet summits
                          (each contains one, or the summit of the max density when absent).


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
PEAK_TYPE="default_peaks"
VAR_WIDTH_PEAKS="0"
WRITE_PVALS=""
SMOOTHING_PARAM=""

# Note: Options in the string that are not immediately followed by ':'
# will not get a value read for them.  Examples are h and P and s.
while getopts 'hc:C:M:e:f:F:m:n:p:S:w:P' opt; do
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
    p)
      PEAK_TYPE=$OPTARG
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
if [ "$PEAK_TYPE" == "default_peaks" ]; then
    require_exes modwt bedGraphToBigWig bedmap samtools hotspot2_part1 hotspot2_part2
else
    if [ "$PEAK_TYPE" == "always_summit_centered" ]; then
	require_exes modwt bedGraphToBigWig bedmap samtools hotspot2_part1 hotspot2_part2 resolveOverlapsInSummit-CenteredPeaks
    else
	require_exes modwt bedGraphToBigWig bedmap samtools hotspot2_part1 hotspot2_part2 resolveOverlapsInSummit-CenteredPeaks findVarWidthPeaks
    fi
fi

CUTCOUNT_EXE="$(dirname "$0")/cutcounts.bash"
DENSPK_EXE="$(dirname "$0")/density-peaks.bash"
MERGE_EXE="$(dirname "$0")/hsmerge.sh"
HOTSPOT_EXE1=hotspot2_part1
HOTSPOT_EXE2=hotspot2_part2

# If an optional peak definition was supplied,
# ensure it's valid.
if [ "$PEAK_TYPE" != "default_peaks" ]; then
    if [ "$PEAK_TYPE" != "always_summit_centered" ]; then
	# value must be varWidth_nn, where nn is an integer
	type=`echo $PEAK_TYPE | cut -f1 -d '_'`
	if [ "$type" != "varWidth" ]; then
	    echo -e "Error in \"-p $PEAK_TYPE\":  \"$PEAK_TYPE\" is invalid."
	    usage
	fi
	minWidthStr=`echo $PEAK_TYPE | cut -f2 -d '_'`
	integerStr=`echo $minWidthStr | awk '{str=sprintf("%d",$1);print str}'`
	if [ "$minWidthStr" == "" ] || [ "$minWidthStr" != "$integerStr" ]; then
	    echo -e "Error in \"-p $PEAK_TYPE\":  \"$PEAK_TYPE\" is invalid."
	    usage
	fi
	id=`echo $PEAK_TYPE | cut -f3- -d '_'`
	if [ "$id" == "" ]; then
	    echo -e "Error in \"-p $PEAK_TYPE\":  \"$PEAK_TYPE\" is invalid."
	    usage
	fi
	VAR_WIDTH_PEAKS="1"
    fi
fi

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
bash "$CUTCOUNT_EXE" "$BAM" "$CUTCOUNTS" "$FRAGMENTS_OUTFILE" "$TOTALCUTS_OUTFILE" "$CHROM_SIZES" $MAPPABLE_REGIONS

if [ ! -s $OUTFILE ] && ([ ! -s $TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1 ] || [ ! -s $TEMP_PVALS ] || [ ! -s $TEMP_CHROM_MAPPING_HOTSPOT2PART1 ]) then
    log "Tallying filtered cut counts in small windows and running part 1 of hotspot2..."
    bedmap --faster --range "$SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE" --delim "\t" --prec 0 --echo --sum "$CENTER_SITES" "$CUTCOUNTS" \
	| "$AWK_EXE" 'BEGIN{OFS="\t"}{if("NAN"==$4){$4=0} print $1, $2, $3, "i", $4}' \
	| "$HOTSPOT_EXE1" --background_size="$BACKGROUND_WINDOW_SIZE" -c $TEMP_CHROM_MAPPING_HOTSPOT2PART1 -p $TEMP_PVALS $SMOOTHING_PARAM -o $TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1
    if [ "$?" != "0" ]; then
	echo -e "An error occurred when calling bedmap on the \"center sites\" and filtered cut counts file, or while running part 1 of hotspot2."
	exit 2
    fi
fi

if [ ! -s $OUTFILE ]; then
    log "Running part 2 of hotspot2..."
    numEntries=`wc -l < $TEMP_PVALS` # used to aid memory allocation
    "$HOTSPOT_EXE2" --fdr_threshold="$CALL_THRESHOLD" $WRITE_PVALS \
       -i $TEMP_INTERMEDIATE_FILE_HOTSPOT2PART1 -n $numEntries -c $TEMP_CHROM_MAPPING_HOTSPOT2PART1 -p $TEMP_PVALS \
       | starch - \
       >"$OUTFILE"
fi

if [ ! -s $HOTSPOT_OUTFILE ]; then
    log "Calling hotspots..."

    "$MERGE_EXE" \
	-f "$HOTSPOT_FDR_THRESHOLD" \
	-m "$MIN_HOTSPOT_WIDTH" \
	"$OUTFILE" \
	"$HOTSPOT_OUTFILE"
fi

if [ ! -s $SPOT_SCORE_OUTFILE ]; then
    log "Calculating SPOT score..."
    num_cleaves=$(cat "$TOTALCUTS_OUTFILE")
    cleaves_in_hotspots=$(bedops --ec -e 1 "$CUTCOUNTS" "$HOTSPOT_OUTFILE" | awk 'BEGIN{s=0} {s+=$5} END {print s}')
    echo "scale=4; $cleaves_in_hotspots / $num_cleaves" \
	| bc \
	>"$SPOT_SCORE_OUTFILE"
fi

if [ ! -s $DENSITY_OUTFILE ] || [ ! -s $PEAKS_OUTFILE ]; then
    #log "Creating peaks and density..."
    if [ "$VAR_WIDTH_PEAKS" == "1" ]; then
	if [ ! -s $TOTALCUTS_OUTFILE ]; then
	    echo -e "Error:  $TOTALCUTS_OUTFILE, which is needed for the detection of variable-width peaks, was not found, or it is empty."
	    exit 2
	fi
	bash "$DENSPK_EXE" "$TMPDIR" "$PEAK_TYPE" "$CUTCOUNTS" "$HOTSPOT_OUTFILE" "$CHROM_SIZES" "$DENSITY_OUTFILE" "$PEAKS_OUTFILE" `cat $TOTALCUTS_OUTFILE`
    else
	bash "$DENSPK_EXE" "$TMPDIR" "$PEAK_TYPE" "$CUTCOUNTS" "$HOTSPOT_OUTFILE" "$CHROM_SIZES" "$DENSITY_OUTFILE" "$PEAKS_OUTFILE"
    fi
fi

if [ ! -s $DENSITY_BW ]; then
    log "Converting density to bigwig..."
    TMPFRAGS="$(mktemp -t fragsXXXXX)"
    unstarch "$DENSITY_OUTFILE" | cut -f1,2,3,5 >"$TMPFRAGS"
    bedGraphToBigWig \
	"$TMPFRAGS" \
	<(cut -f1,3 "$CHROM_SIZES") \
	"$DENSITY_BW"
fi

log "Done!"

rm -f "$TMPFRAGS"
if [[ $clean != 0 ]]; then
  rm -rf "$TMPDIR"
fi

exit 0
