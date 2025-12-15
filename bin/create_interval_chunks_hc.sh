#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = counts_per_chunk
# $4 = split_overweight
# $5 = counts_file

TARGET_COUNTS_PER_CHUNK=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')
OUTDIR=$(pwd)
SPLIT_OVERWEIGHT=${4}

# HC interval chunks operate on a single sample counts file only 

# Ensure exactly 4 columns (chr, start, end, count) and drop count==0
awk 'NF>=4 && $4+0 != 0 {print $1"\t"$2"\t"$3"\t"$4}' "${5}" > intervals_with_counts.bed

# Optionally Split intervals that individually exceed the target counts.
# Assumes counts are roughly uniform across the interval length.
# WARNING: Makes more even intervals at risk of artefacts near interval end.

if [[ "$SPLIT_OVERWEIGHT"  == "true" ]]; then
  awk -v target="$TARGET_COUNTS_PER_CHUNK" '
  BEGIN{OFS="\t"}
  {
    chrom=$1; start=$2; end=$3; w=$4
    len=end-start

    # guard against weird zero/negative lengths
    if (len <= 0) next

    # if not overweight, keep as-is
    if (w <= target) {
      print chrom, start, end, w
      next
    }

    # number of pieces needed (ceil)
    n = int((w + target - 1) / target)

    # don’t create more pieces than bases
    if (n > len) n = len

    base = int(len / n)
    rem  = len - base * n

    substart = start
    for (i=1; i<=n; i++) {
      sz = base + (i <= rem ? 1 : 0)
      subend = substart + sz

      # proportional weight by length to preserve density
      subw = w * sz / len

      # round to integer
      subw = int(subw + 0.5)
      print chrom, substart, subend, subw
      substart = subend
    }
  }
  ' intervals_with_counts.bed > intervals_split.bed
else
  cat intervals_with_counts.bed > intervals_split.bed
fi

# Calculate total counts and number of intervals
TOTAL_COUNTS=$( awk '{s+=$4} END{print s+0}' intervals_split.bed )
N_INTERVALS=$(cat intervals_split.bed | wc -l)

# Decide number of file splits (chunks) to keep counts balanced.
K=$(( (TOTAL_COUNTS + TARGET_COUNTS_PER_CHUNK - 1) / TARGET_COUNTS_PER_CHUNK )) 
if [ "$K" -lt 1 ]; then K=1; fi
if [ "$K" -gt "$N_INTERVALS" ]; then K="$N_INTERVALS"; fi

# Assign intervals to chunks, keep it contiguous so they are in sorted order
awk -v outdir="$OUTDIR" -v tot="$TOTAL_COUNTS" -v n="$N_INTERVALS" -v K="$K" '
function abs(x){ return x<0 ? -x : x }

BEGIN{
  chunk=1
  remK=K
  remW=tot
  sum=0
  ideal = remW / remK
  fname = sprintf("%s/chunk_%d.bed", outdir, chunk)
}

{
  chrom=$1; start=$2; end=$3; w=$4+0
  lines_after = n - FNR

  if (sum>0 && remK>1) {
    stop_now_better = (abs(sum-ideal) <= abs((sum+w)-ideal))
    enough_left     = (lines_after >= (remK-2))

    if (stop_now_better && enough_left) {
      chunk++
      remK--
      sum=0
      ideal = remW / remK
      fname = sprintf("%s/chunk_%d.bed", outdir, chunk)
    }
  }

  print chrom"\t"start"\t"end"\t"w > fname
  sum  += w
  remW -= w
}
' intervals_split.bed
   
# Rename each output file
for i in *chunk_*.bed;do
  # Pad output chunk names
  n=$(basename "$i" | sed -E 's/^chunk_([0-9]+)\.bed/\1/')
  pad=$(printf "%05d" "$n")

  # compute hash of this chunk’s contents
  hash=$(md5sum "$i" | awk '{print $1}')

  # merge any abutting intervals
  out="_${pad}${hash}.bed"
  cut -f1-4 "$i" \
  | bedtools merge -i stdin > "$out"
  
  # report size
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$i")
  counts=$(awk '{s+=$4} END{print s+0}' "$i")
  echo "${out}: ${bases} genomic bases, ${counts} summed counts"

  rm $i
done
