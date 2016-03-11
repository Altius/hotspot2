#!/bin/bash
set -x -e -o pipefail

if [[ $# != 5 ]] ; then
  echo "Usage: $0 tmpdir wavelets-binary chrom-sizes.bed density-out.starch peaks-out.starch" >&2
  exit 2
fi

tmpdir=$1
wavelets=$2
chrfile=$3
density=$4
pk=$5

if [ ! -d $tmpdir ] then
  mkdir -p $tmpdir
fi

## density params
bins=150
step=20
halfbin=$((bins/2))
rangepad=$((bins/2-step/2))

## wavelet peakfinding params
wavletlvl=3
filter=Haar
boundary_type=reflected
lvl=3


## Tag density, 150bp window, sliding every 20bp, used for peak-finding.
echo "calculating densities..."
sort-bed $chrfile \
  | awk -v b=$bins -v s=$step \
     'BEGIN {OFS="\t"; hs=s/2; hb=b/2} ; { \
       for ( start = $2+hb-hs; start < $3-hb-hs; start+=s) { \
         print $1, start, start+s, "."; \
       } \
     }' \
  | bedmap --faster --range $rangepad --delim "\t" --echo --count - $tagsb \
  | starch - \
 > $density

echo "peak-finding..."
outs=""
for chr in $(sort-bed $chrfile | cut -f1)
do
  unstarch $chr $density \
    | cut -f5 \
    | wavelets --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter - \
    > $tmpdir/.waves

  unstarch $chr $density \
    | paste - $tmpdir/.waves \
    | cut -f1-4,6 \
    | awk 'BEGIN{incr=0} ; { \
            if ( NR > 0 ) { \
              if ( incr == 1 ) { \
                if ( $5-lastv < 0 ) { \
                  print lastl; incr=0; \
                } \
              } else { \
                if ( $5-lastv > 0 ) incr=1; \
              } \
            } \
            lastv=$5; lastl=$0; \
          }' \
    > $tmpdir/.wave.$chr

  outs="$outs $tmpdir/.wave.$chr"
done

sort-bed $outs \
  | awk -v h=$halfbin '{m=($2+$3)/2; left=m-h; if(left < 0) left=0; print $1"\t"left"\t"m+h"\t"$4"\t"$5}' - \
  > $pk

exit 0