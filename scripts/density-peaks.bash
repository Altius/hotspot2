#!/bin/bash
set -x -e -o pipefail

log() {
  echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t$*"
}

if [[ $# != 7 ]]; then
  echo "Usage: $0 tmpdir wavelets-binary <tags.starch> <hotspots.starch> <chrom-sizes.bed> <density-out.starch> <peaks-out.starch>" >&2
  exit 2
fi

tmpdir=$1
wavelets=$2
tags=$3
hotspots=$4
chrfile=$5
density=$6
pk=$7

if [ ! -d "$tmpdir" ]; then
  mkdir -p "$tmpdir"
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

log "Calculating densities and peak-finding..."
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
    | $wavelets --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter_type - \
     >"$tmpdir/.waves"

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

  unstarch "$chr" "$tmpdir/.dens.$chr.starch" \
      | paste - "$tmpdir/.waves" \
      | "$AWK_EXE" 'BEGIN{OFS="\t"; increasing=0} \
            {if (NR > 0) { \
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
            }END{if(increasing && NR>0){split(prevLine, x, "\t"); print x[1],x[2],x[3],x[4],x[5]}}' \
      | bedmap --ec --sweep-all --chrom "$chr" --skip-unmapped --fraction-map 1.0 --echo --echo-map $hotspots - \
      | "$AWK_EXE" -F "|" 'BEGIN{OFS="\t"}{split($1,y,"\t");lengthOfX=split($2,x,";"); \
                                  for(i=1;i<=lengthOfX;i++){ \
                                     split(x[i],z,"\t"); \
                                     print y[1],y[2],y[3],z[5],int((z[2]+z[3])/2),int((y[2]+y[3])/2) \
                                  } \
                                }' \
      | closest-features --chrom "$chr" --no-overlaps - $hotspots \
      | "$AWK_EXE" -v chromLength=$chrLength -v "halfWidth=$halfbin" -F "|" \
           'BEGIN{OFS = "\t"} \
            {if("NA" == $2) { \
                if(1 == NR) { L_xR = 0 } \
                else { L_xR = prevL_xR } \
             } \
             else { \
                split($2,z,"\t"); \
                if(NR > 0 && prevL_xR > z[3]) { \
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
      > "$tmpdir/.wave-pks.$bins.inHotspots.$chr"

  # We expect, and want, every hotspot to contain at least one peak.
  # Some hotspots won't have any wavelet peaks called within them;
  # for these hotspots, we force a peak call within each one,
  # using the 20bp within it that contain the maximum density value
  # across the hotspot as a starting point.

  # When a hotspot is sufficiently wide, pad the 20-bp "max region" symmetrically.
  # Otherwise, if the hotspot is wide enough to fully contain a peak, pad the 20-bp "max element" asymmetrically,
  # with its left (or right) boundary at the left (or right) boundary of the enclosing hotspot.
  # And if the hotspot is not wide enough to fully contain a peak,
  # if the region between the right boundary of the upstream hotspot (L_xR)
  # and the left boundary of the downstream hotspot (R_xL) is wide enough to fully contain a peak,
  # place it symmetrically across the hotspot, or asymmetrically (boundary at L_xR or R_xL) when necessary.
  # And if the region between L_xR and R_xL is too narrow to fully contain a peak,
  # define the peak to begin at L_xR and end at R_xL.
  # To help flag forced peak calls as such,
  # we report the score of each with ".0" after it (i.e., with precision = 1).

  bedops --ec --chrom $chr -n 1 $hotspots "$tmpdir/.wave-pks.$bins.inHotspots.$chr" \
      | bedmap --prec 1 --echo --max-element - "$tmpdir/.dens.$chr.starch" \
      | "$AWK_EXE" -F "|" 'BEGIN{OFS="\t"}{split($1,y,"\t");split($2,z,"\t"); \
                                     print y[1],y[2],y[3],z[5],int((z[2]+z[3])/2),int((y[2]+y[3])/2); \
                                }' \
      | closest-features --chrom "$chr" --no-overlaps - "$tmpdir/.wave-pks.$bins.inHotspots.$chr" \
      | "$AWK_EXE" -v chromLength=$chrLength -v "halfWidth=$halfbin" -F "|" \
           'BEGIN{OFS = "\t"} \
            {if("NA" == $2) { \
                L_xR = 0; \
             } \
             else { \
                split($2,z,"\t"); \
                if(NR > 0 && prevL_xR > z[3]) { \
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
      > "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"

  namesOfNonemptyPeakFiles=""
  if [ -s "$tmpdir/.wave-pks.$bins.inHotspots.$chr" ]; then
     namesOfNonemptyPeakFiles="$tmpdir/.wave-pks.$bins.inHotspots.$chr"
  fi

  if [ -s "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr" ]; then
     namesOfNonemptyPeakFiles="$namesOfNonemptyPeakFiles $tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr"
  fi
  if [ "" != "$namesOfNonemptyPeakFiles" ]; then
     bedops -u $namesOfNonemptyPeakFiles \
      >"$tmpdir/.wave-pks.$bins.finalPeaks.$chr"
     pkouts="$pkouts $tmpdir/.wave-pks.$bins.finalPeaks.$chr"
  fi

  rm -f "$tmpdir/.wave-pks.$bins.inHotspots_uniq.$chr" "$tmpdir/.wave-pks.$bins.inHotspots_newMerged.$chr" "$tmpdir/.wave-pks.$bins.forcedPeakCalls.$chr" \

done

if [ "" == "$pkouts" ]; then
  log "Finishing; no peaks were called(!)..."
else
  log "Finalizing peaks..."
  cat $pkouts \
    | bedmap --echo --skip-unmapped --sweep-all --fraction-either 0.25 - "$hotspots" \
    | starch - \
      >"$pk"

  unstarch "$pk" \
    | awk -v "h=$halfbin" 'BEGIN {OFS="\t"} ; { print $1, $2, $3, ".", "0", ".", $5, "-1", "-1", h }' \
    | starch - \
      >"${pk/.starch/.narrowpeaks.starch}"
fi

log "Finalizing density..."
starchcat $densouts \
  >"$density"

exit 0
