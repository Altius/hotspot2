#!/bin/bash
set -x -e -o pipefail

log() {
  echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t$*"
}

usage () {
  echo "Usage: $0 tmpdir peak_definition <tags.starch> <hotspots.starch> <chrom-sizes.bed> <density-out.starch> <peaks-out.starch> [tally of mapped cleavages]" >&2
  echo "where peak_definition is either \"default_peaks\", \"always_summit_centered\", or \"varWidth_nn_ID\"," >&2
  echo "with these exact capitalizations (or lack thereof)." >&2
  echo "In \"varWidth_nn_ID\", nn is an integer specifying the minimum width (in bp) required for a variable-width peak to be reported," >&2
  echo "and ID is an identifier that will be written into column 4 of every line of the peaks output file." >&2
  echo "If \"varWidth_nn_ID\" is received, then the total number of mapped cleavages must also be supplied as the final argument to $0." >&2
  exit 2
}

if [[ $# != 7 ]] && [[ $# != 8 ]]; then
   usage
fi

tmpdir=$1
peak_type=$2
tags=$3
hotspots=$4
chrfile=$5
density=$6
pk=$7

if [ ! -d "$tmpdir" ]; then
  mkdir -p "$tmpdir"
fi

# Parse the peak definition.
VAR_WIDTH_PEAKS="0"
MIN_WIDTH=""
VAR_WIDTH_PEAKS_ID=""
if [ "$peak_type" != "default_peaks" ]; then
    if [ "$peak_type" != "always_summit_centered" ]; then
	# value must be varWidth_nn, where nn is an integer
	type=`echo $peak_type | cut -f1 -d '_'`
	if [ "$type" != "varWidth" ]; then
	    echo -e "Error in \"-p $peak_type\":  \"$peak_type\" is invalid."
	    usage
	fi
	minWidthStr=`echo $peak_type | cut -f2 -d '_'`
	integerStr=`echo $minWidthStr | awk '{str=sprintf("%d",$1);print str}'`
	if [ "$minWidthStr" == "" ] || [ "$minWidthStr" != "$integerStr" ]; then
	    echo -e "Error in \"-p $peak_type\":  \"$peak_type\" is invalid."
	    usage
	fi
	MIN_WIDTH=$minWidthStr
	VAR_WIDTH_PEAKS_ID=`echo $peak_type | cut -f3- -d '_'`
	if [ "$8" != "" ]; then
	    TOTAL_TAGS=$8
	else
	    echo -e "Error:  Peak definition \"$peak_type\" was received by $0,"
	    echo -e "but required accompanying 8th argument (total mapped cleavages) was not supplied."
	    usage
	fi
	VAR_WIDTH_PEAKS="1"
    fi
fi

# Get the required executables.
# The script that called this script (hotspot2.sh) already ensured that `which modwt` succeeded.
# We double-check these here just in case the user called the current script directly.
get_wavelets=`which modwt 2> /dev/null`
if [ ! -x $get_wavelets ]; then
   echo -e "Error:  Required executable \"modwt\" was not found, or permission to execute it was not found."
   exit 2
fi
PEAK_OVERLAP_RESOLVER=""
FIND_VAR_WIDTH_PEAKS_EXE=""
if [ "$peak_type" != "default_peaks" ]; then
    PEAK_OVERLAP_RESOLVER=`which resolveOverlapsInSummit-CenteredPeaks 2> /dev/null`
    if [ ! -x $PEAK_OVERLAP_RESOLVER ]; then
	echo -e "Error:  Required executable \"resolveOverlapsInSummit-CenteredPeaks\" was not found, or permission to execute it was not found."
	exit 2
    fi
    if [ "$VAR_WIDTH_PEAKS" == "1" ]; then
	FIND_VAR_WIDTH_PEAKS_EXE=`which findVarWidthPeaks 2> /dev/null`
	if [ ! -x $FIND_VAR_WIDTH_PEAKS_EXE ]; then
	    echo -e "Error:  Required executable \"findVarWidthPeaks\" was not found, or permission to execute it was not found."
	    exit 2
	fi	
    fi
fi

## density params
bins=150
step=20
halfbin=$((bins / 2))
rangepad=$((bins / 2 - step / 2))

## wavelet peakfinding params
waveletlvl=3
filter_type=Haar
boundary_type=reflected

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)

# Function to compute wavelet peaks in the "standard" case,
# in which each reported peak is expected to lie entirely within a hotspot.
# (Exception:  Hotspots narrower than the (fixed) width of the peak.)
#------------------------------------------
  # Get padded wavelet peaks that overlap hotspots.
  # Write the density value as the "score" for each wavelet peak.
  # Require all of each 20-bp wavelet peak to land within a hotspot to get called as a peak for that hotspot.
  # When a hotspot is sufficiently wide, pad the 20-bp wavelet peak symmetrically.
  # Otherwise, if the hotspot is wide enough to fully contain a peak, pad the 20-bp wavelet peak asymmetrically,
  # with its left (or right) boundary at the left (or right) boundary of the enclosing hotspot.
  # And if the hotspot is not wide enough to fully contain a peak,
  # if the region between the right boundary of the upstream hotspot (L_xR)
  # and the left boundary of the downstream hotspot (R_xL) is wide enough to fully contain a peak,
  # place it symmetrically across the hotspot, or asymmetrically (boundary at L_xR or R_xL) when necessary.
  # And if the region between L_xR and R_xL is too narrow to fully contain a peak,
  # define the peak to begin at L_xR and end at R_xL.
#
# Arguments:
# $1 = chrom
# $2 = chromLength
# $3 = densityFile
# $4 = .waves file
# $5 = hotspots
function get_wavelet_peaks_standard {
    unstarch $1 $3 \
	| paste - $4 \
        | "$AWK_EXE" 'BEGIN{OFS="\t"; increasing=0} \
            {if (NR > 1) { \
                if (increasing) { \
                   if ($6 < prevValue) { \
                      split(prevLine, x, "\t"); \
                      print x[1], x[2], x[3], x[4], x[5]; \
                      increasing = 0; \
                   } \
                } \
                else { \
                   if ($6 > prevValue) increasing = 1; \
               } \
             } \
             prevValue=$6; prevLine=$0; \
            }END{if(increasing){split(prevLine, x, "\t"); print x[1],x[2],x[3],x[4],x[5]}}' \
     | bedmap --ec --sweep-all --chrom $1 --skip-unmapped --fraction-map 1.0 --echo --echo-map $5 - \
     | "$AWK_EXE" -F "|" 'BEGIN{OFS="\t"}{split($1,y,"\t");lengthOfX=split($2,x,";"); \
                                  for(i=1;i<=lengthOfX;i++){ \
                                     split(x[i],z,"\t"); \
                                     print y[1],y[2],y[3],z[5],int((z[2]+z[3])/2),int((y[2]+y[3])/2) \
                                  } \
                                }' \
     | closest-features --chrom $1 --no-overlaps - $5 \
     | "$AWK_EXE" -v chromLength="$2" -v "halfWidth=$halfbin" -F "|" \
           'BEGIN{OFS = "\t"} \
            {if("NA" == $2) { \
                if(1 == NR) { L_xR = 0 } \
                else { L_xR = prevL_xR } \
             } \
             else { \
                split($2,z,"\t"); \
                if(NR > 1 && prevL_xR > z[3]) { \
                   L_xR = prevL_xR; \
                } \
                else { \
                   L_xR = z[3]; \
                } \
             } \
             if("NA" == $3) { \
                R_xL = chromLength; \
             } \
             else { \
                split($3,z,"\t"); \
                R_xL = z[2]; \
             } \
             split($1, y, "\t"); \
             chr = y[1]; \
             C_xL = y[2]; \
             C_xR = y[3]; \
             score = y[4]; \
             x_centerOfWavelet = y[5]; \
             x_centerOfHotspot = y[6]; \
             if (L_xR > C_xL) { C_xL = L_xR } \
             if (C_xR - C_xL >= 2*halfWidth) { \
                if (C_xL > x_centerOfWavelet - halfWidth) { \
                   xL = C_xL; \
                   xR = xL + 2*halfWidth; \
                } \
                else { \
                   if (C_xR < x_centerOfWavelet + halfWidth) { \
                      xR = C_xR; \
                      xL = xR - 2*halfWidth; \
                   } \
                   else { \
                      xL = x_centerOfWavelet - halfWidth; \
                      xR = x_centerOfWavelet + halfWidth; \
                   } \
                } \
             } \
             else { \
                if (R_xL - L_xR >= 2*halfWidth) { \
                   if (L_xR > x_centerOfHotspot - halfWidth) { \
                      xL = L_xR; \
                      xR = L_xR + 2*halfWidth; \
                   } \
                   else { \
                      if (R_xL < x_centerOfHotspot + halfWidth) { \
                         xR = R_xL; \
                         xL = xR - 2*halfWidth; \
                      } \
                      else { \
                         xL = x_centerOfHotspot - halfWidth; \
                         xR = x_centerOfHotspot + halfWidth; \
                      } \
                   } \
                } \
                else { \
                   xL = L_xR; \
                   xR = R_xL; \
                   if (xL >= xR) { \
                      prevL_xR = L_xR; \
                      next; \
                   } \
                } \
             } \
             print chr, xL, xR, "i", score; \
             prevL_xR = xR; \
            }' \
     > $6
}
#------------------------------------------


# Function to compute wavelet peaks in the "optional" case,
# in which each reported 150-bp peak is centered on the central bp of the 20-bp summit bin
# (i.e., no shifting to attempt to make the 150-bp peak lie entirely within a hotspot,
# which is performed in the "standard" case).
#------------------------------------------
#
# Arguments:
# $1 = chrom
# $2 = densityFile
# $3 = .waves file
# $4 = hotspots
# $5 = outfile
function get_wavelet_peaks_always_centered_on_summits {
  unstarch $1 $2 \
      | paste - $3 \
      | "$AWK_EXE" -v halfWidth=$halfbin 'BEGIN{OFS="\t"; increasing=0} \
            {if (NR > 1) { \
                if (increasing) { \
                   if ($6 < prevValue) { \
                      split(prevLine, x, "\t"); \
                      print x[1], x[2], x[3], "i", x[5]; \
                      increasing = 0; \
                   } \
                } \
                else { \
                   if ($6 > prevValue) increasing = 1; \
                } \
             } \
             prevValue=$6; prevLine=$0; \
            }END{if(increasing){split(prevLine, x, "\t");print x[1], x[2], x[3], "i", x[5]}}' \
      | bedmap --faster --chrom $1 --skip-unmapped --echo --echo-map - $4 \
      | "$AWK_EXE" -F "|" -v halfWidth=$halfbin 'BEGIN{OFS="\t"}{split($1,x,"\t"); \
            chr=x[1]; score=x[5]; xC=int(0.5*(x[2]+x[3])); \
            numOverlappingHotspots = split($2,y,";"); \
            if (2 == numOverlappingHotspots) { \
               split(y[1],x,"\t"); L_xR=x[3]; \
               split(y[2],z,"\t"); R_xL=z[2]+1; \
               if (xC >= R_xL || xC <= L_xR) { \
                  print chr, xC-halfWidth, xC+halfWidth, "i", score; \
               } \
               else { \
                  distL = xC - L_xR; \
                  distR = R_xL - xC; \
                  if (distR < distL) { \
                     print chr, R_xL-halfWidth, R_xL+halfWidth, "i", score; \
                  } \
                  else { \
                     print chr, L_xR-halfWidth, L_xR+halfWidth, "i", score; \
                  } \
               } \
            } \
            else { \
               split(y[1],x,"\t"); hsBeg=x[2]+1; hsEnd=x[3]; \
               if (xC < hsBeg) { \
                  print chr, hsBeg-halfWidth, hsBeg+halfWidth, "i", score; \
               } \
               else { \
                  if (xC <= hsEnd) { \
                     print chr, xC-halfWidth, xC+halfWidth, "i", score; \
                  } \
                  else { \
                     print chr, hsEnd-halfWidth, hsEnd+halfWidth, "i", score; \
                  } \
               } \
            } \
           }' \
      > $5
}
#------------------------------------------


# Function to force-call peaks (one per hotspot for which no wavelet peak was found)
# in the "standard" case.
#------------------------------------------
  # We expect, and want, every hotspot to contain at least one peak.
  # Some hotspots won't have any wavelet peaks called within them;
  # for these hotspots, we force a peak call within each one,
  # using the 20bp within it that contain the maximum density value
  # across the hotspot as a starting point.
# Arguments:
# $1 = chrom
# $2 = chromLength
# $3 = wavePeaksInHotspots
# $4 = densityFile
# $5 = hotspots
# $6 = outfile
function get_forcedCall_peaks_standard {
  if [ -s $3 ]; then
     bedops --ec --chrom $1 -n 1 $5 $3 \
         | bedmap --prec 1 --echo --max-element - $4 \
         | "$AWK_EXE" -F "|" 'BEGIN{OFS="\t"}{split($1,y,"\t");split($2,z,"\t"); \
                                    print y[1],y[2],y[3],z[5],int((z[2]+z[3])/2),int((y[2]+y[3])/2); \
                              }' \
         | closest-features --chrom $1 --no-overlaps - $3 \
         | "$AWK_EXE" -v chromLength=$2 -v "halfWidth=$halfbin" -F "|" \
           'BEGIN{OFS = "\t"} \
            {if("NA" == $2) { \
                L_xR = 0; \
             } \
             else { \
                split($2,z,"\t"); \
                if(NR > 1 && prevL_xR > z[3]) { \
                   L_xR = prevL_xR; \
                } \
                else { \
                   L_xR = z[3]; \
                } \
             } \
             if("NA" == $3) { \
                R_xL = chromLength; \
             } \
             else { \
                split($3,z,"\t"); \
                R_xL = z[2]; \
             } \
             split($1, y, "\t"); \
             chr = y[1]; \
             C_xL = y[2]; \
             C_xR = y[3]; \
             score = y[4]; \
             x_maxSignal = y[5]; \
             x_centerOfHotspot = y[6]; \
             if (L_xR > C_xL) { C_xL = L_xR } \
             if (L_xR > C_xR) { \
                # this entire hotspot is fully enclosed by the previous peak call! \
                prevL_xR = L_xR; \
                next; \
             } \
             if (C_xR - C_xL >= 2*halfWidth) { \
                if (C_xL > x_maxSignal - halfWidth) { \
                   xL = C_xL; \
                   xR = xL + 2*halfWidth; \
                } \
                else { \
                   if (C_xR < x_maxSignal + halfWidth) { \
                      xR = C_xR; \
                      xL = xR - 2*halfWidth; \
                   } \
                   else { \
                      xL = x_maxSignal - halfWidth; \
                      xR = x_maxSignal + halfWidth; \
                   } \
                } \
             } \
             else { \
                if (R_xL - L_xR >= 2*halfWidth) { \
                   if (L_xR > x_centerOfHotspot - halfWidth) { \
                      xL = L_xR; \
                      xR = L_xR + 2*halfWidth; \
                   } \
                   else { \
                      if (R_xL < x_centerOfHotspot + halfWidth) { \
                         xR = R_xL; \
                         xL = xR - 2*halfWidth; \
                      } \
                      else { \
                         xL = x_centerOfHotspot - halfWidth; \
                         xR = x_centerOfHotspot + halfWidth; \
                      } \
                   } \
                } \
                else { \
                   xL = L_xR; \
                   xR = R_xL; \
                } \
                if (xR > C_xR && C_xR-xL < 0.25*halfWidth && L_xR > C_xL) { \
                   # if only an eighth or less of this peak actually lies within this hotspot, \
                   # and the previous peak (for the previous hotspot) \
                   # overlaps the current hotspot by at least 1bp, \
                   # do not bother to call another peak for this hotspot \
                   prevL_xR = L_xR; \
                   next; \
                } \
             } \
             print chr, xL, xR, "i", score; \
             prevL_xR = xR; \
            }' \
         > $6
  fi
}
#------------------------------------------


# Function to force-call peaks (one per hotspot for which no wavelet peak was found)
# in the optional "always centered on the summit" case.
#------------------------------------------
# Arguments:
# $1 = chrom
# $2 = wavePeaksInHotspots
# $3 = densityFile
# $4 = hotspots
# $5 = outfile
function get_forcedCall_peaks_always_centered_on_summits {
  if [ -s $2 ]; then
     bedops --ec --chrom $1 -n 1 $4 $2 \
         | bedmap --faster --chrom $1 --prec 1 --echo --max-element - $3 \
         | "$AWK_EXE" -F "|" -v halfWidth=$halfbin 'BEGIN{OFS="\t"}{split($1,x,"\t"); \
              chr=x[1]; hsBeg=x[2]+1; hsEnd=x[3]; \
              split($2,x,"\t"); xC = int(0.5*(x[2]+x[3])); score=x[5]; \
              if (xC < hsBeg) { \
                 print chr, hsBeg-halfWidth, hsBeg+halfWidth, "i", score; \
              } \
              else { \
                 if (xC <= hsEnd) { \
                    print chr, xC-halfWidth, xC+halfWidth, "i", score; \
                 } \
                 else { \
                    print chr, hsEnd-halfWidth, hsEnd+halfWidth, "i", score; \
                 } \
              } \
             }' \
	 > $5
  fi
}
#------------------------------------------


log "Calculating densities and finding peaks..."
pkouts=""
densouts=""
for chr in $(awk '{print $1}' "$chrfile"); do
  log "\tProcessing $chr"

  chrLength=`grep -w ^"$chr" "$chrfile" | cut -f3`

  ## Tag density, 150bp window, sliding every 20bp, used for peak-finding and display.
  ##  --sweep-all used to prevent a possible broken pipe.
  ## Use awk, not sed to change NAN to 0 since a 'chromosome' name could include NAN and then fields 1 and 4 would change.
  ## sed could be used with regular expressions, but awk is easy enough here, I think.

  bedops --ec --chop 20 --stagger 20 --chrom "$chr" "$chrfile" \
    | bedmap --faster --sweep-all --chrom "$chr" --range "$rangepad" --delim "\t" --prec 0 --echo --echo-ref-row-id --sum - "$tags" \
    | awk 'BEGIN {OFS="\t"} ; { if ( $(NF) == "NAN" ) { $(NF)=0; } print; }' \
    | starch - \
     >"$tmpdir/.dens.$chr.starch"

  densouts="$densouts $tmpdir/.dens.$chr.starch"

  unstarch "$chr" "$tmpdir/.dens.$chr.starch" \
    | cut -f5 \
    | $get_wavelets --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter_type - \
     >"$tmpdir/.waves"

  # Get the wavelet peaks for this chromosome.  (Function definitions are given earlier in this script.)
  if [ "$peak_type" == "default_peaks" ]; then
     get_wavelet_peaks_standard $chr $chrLength "$tmpdir/.dens.$chr.starch" "$tmpdir/.waves" $hotspots "$tmpdir/.wave-pks.$bins.inHotspots.$chr"
     if [ $? != "0" ]; then
        echo "An error occurred during get_wavelet_peaks_standard() while processing $chr in $0."
	exit 2
     fi
  else
     get_wavelet_peaks_always_centered_on_summits $chr "$tmpdir/.dens.$chr.starch" "$tmpdir/.waves" $hotspots "$tmpdir/.wave-pks.$bins.inHotspots.$chr"
     if [ $? != "0" ]; then
        echo "An error occurred during get_wavelet_peaks_always_centered_on_summits() while processing $chr in $0."
	exit 2
     fi
  fi
  
  # Get "forced-call" peaks for this chromosome.  (Explanation and function definitions are given earlier in this script.)
  if [ "$peak_type" == "default_peaks" ]; then
     get_forcedCall_peaks_standard $chr $chrLength "$tmpdir/.wave-pks.$bins.inHotspots.$chr" \
				   "$tmpdir/.dens.$chr.starch" $hotspots "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"
     if [ $? != "0" ]; then
        echo "An error occurred during get_forcedCall_peaks_standard() while processing $chr in $0."
	exit 2
     fi
  else
     get_forcedCall_peaks_always_centered_on_summits $chr "$tmpdir/.wave-pks.$bins.inHotspots.$chr" "$tmpdir/.dens.$chr.starch" \
				                     $hotspots "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"
     if [ $? != "0" ]; then
        echo "An error occurred during get_forcedCall_peaks_always_centered_on_summits() while processing $chr in $0."
	exit 2
     fi
  fi
  
  namesOfNonemptyPeakFiles=""
  if [ -s "$tmpdir/.wave-pks.$bins.inHotspots.$chr" ]; then
     namesOfNonemptyPeakFiles="$tmpdir/.wave-pks.$bins.inHotspots.$chr"
  fi

  if [ -s "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr" ]; then
     namesOfNonemptyPeakFiles="$namesOfNonemptyPeakFiles $tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"
  fi
  if [ "" != "$namesOfNonemptyPeakFiles" ]; then
     if [ "$peak_type" == "default_peaks" ]; then
        bedops -u $namesOfNonemptyPeakFiles \
            >"$tmpdir/.wave-pks.$bins.finalPeaks.$chr"
     else
	 # In this case, there are wavelet and forced-call peaks that can overlap.
	 # An external program resolves overlaps, retaining the peaks with the highest signal strengths.
	 # Example, four peaks, labeled A-D, with signal strengths 1-4:
         #
         #               ------- D (4)
         #           ------- C (3)
         #     ------- B (2)
         # ------- A (1)
         #
         # The result will be:
         #
         #               ------- D (4)
         #     ------- B (2)
	 #

	 TEMPFILE="$tmpdir/.wave-pks.$bins.finalPeaksWithOverlaps.$chr"
	 bedops -u $namesOfNonemptyPeakFiles \
	    > $TEMPFILE
	 bedops -m $TEMPFILE \
	     | bedmap --echo --echo-map - $TEMPFILE \
	     | $PEAK_OVERLAP_RESOLVER \
	     >"$tmpdir/.wave-pks.$bins.finalPeaks.$chr"
	 rm -f $TEMPFILE
     fi
     pkouts="$pkouts $tmpdir/.wave-pks.$bins.finalPeaks.$chr"
  fi

  rm -f "$tmpdir/.wave-pks.$bins.inHotspots.$chr" "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"

done

log "Finalizing density..."
starchcat $densouts \
  >"$density"

if [ "$VAR_WIDTH_PEAKS" == "1" ]; then
    TEMPFILE="$tmpdir/.tempVarWidthPeaks"
    bedops -u $pkouts \
       | awk 'BEGIN{OFS="\t"}{c=int(0.5*($2+$3));print $1,c-1,c}' \
       | bedmap --echo --echo-map - $hotspots \
       | awk -F "|" '{print $2"|"$1}' \
       | bedmap --echo --echo-map - $density \
       | awk -F "|" -v totalTags=$TOTAL_TAGS -v ID=$VAR_WIDTH_PEAKS_ID \
           'BEGIN{OFS="\t";c=1000000./totalTags} \
           {split($1,x,"\t"); \
            chr=x[1];hsBeg=x[2];hsEnd=x[3]; \
            split($2,z,"\t"); \
            summit=z[3]; \
            len=split($3,y,";"); \
            for(i=1;i<=len;i++){ \
               split(y[i],x,"\t"); \
               for(s=x[2]; s<x[3]; s++){ \
                  if(s >= hsBeg && s+1 <= hsEnd){ \
                     print chr, s, s+1, ID, c*x[5], summit; \
                  } \
               } \
            } \
           }' \
       | $FIND_VAR_WIDTH_PEAKS_EXE $MIN_WIDTH \
       | sort-bed - \
       | starch - \
       > $TEMPFILE
    pkouts=$TEMPFILE
    if [ ! -s $pkouts ]; then
	pkouts=""
    fi
fi

if [ "" == "$pkouts" ]; then
  log "Finishing; no peaks were called(!)..."
else
  log "Finalizing peaks..."
  if [ "$VAR_WIDTH_PEAKS" == "1" ]; then
    mv $pkouts $pk # $pkouts has already been starched
  else
    bedops -u $pkouts \
      | starch - \
      > $pk
  fi
  unstarch "$pk" \
    | awk -v "h=$halfbin" 'BEGIN {OFS="\t"} ; { print $1, $2, $3, ".", "0", ".", $5, "-1", "-1", h }' \
    | starch - \
      >"${pk/.starch/.narrowpeaks.starch}"
fi

exit 0
