#!/usr/bin/env awk

BEGIN {OFS="\t"; chr=""}
{
  if (x == $4 && y == $5 && z == $6 && end == $2 && chr == $1) {
    end = $3;
  } else {
    if (chr != ""){
      print chr, start, end, x, y, z;
    }
    chr = $1; start = $2; end = $3; x = $4; y = $5; z = $6;
  }
}
END { if (chr != "") {print chr, start, end, x, y, z}}

# NR == 1 {OFS="\t"; chr=""}
# NR > 1 {
#   if (x == $4 && end == $2 && chr == $1) {
#     end = $3;
#   } else {
#     if (chr != ""){
#       print chr, start, end, x;
#     }
#     chr = $1; start = $2; end = $3; x = $4;
#   }
# }
# END { if (chr != "") {print chr, start, end, x}}
