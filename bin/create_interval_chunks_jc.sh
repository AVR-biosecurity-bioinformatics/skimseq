#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = counts_per_chunk
# $4 = counts_files

COUNTS_PER_CHUNK=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')
OUTDIR=$(pwd)

# counts files for all samples
bedtools unionbedg -i *counts.bed -filler 0 > combined_counts.bed


# Take the sum of feature counts across windows
awk 'BEGIN{OFS="\t"} {
  sum = 0
  for (i = 4; i <= NF; i++) sum += $i
  print $1, $2, $3, sum
}' combined_counts.bed > intervals_with_counts.bed

# Use greedy algorithm to assign intervals to chunks.
# Once they reach COUNTS_PER_CHUNK, begin a new chunk.
awk -v target="$COUNTS_PER_CHUNK" -v outdir="$OUTDIR" '
    BEGIN{chunk=1; sum=0; fname=sprintf("%s/chunk_%d.bed", outdir, chunk)}
    {
    chrom=$1; start=$2; end=$3; weighted=$4
    if(sum+weighted>target && sum>0){
        chunk++
        fname=sprintf("%s/chunk_%d.bed", outdir, chunk)
        sum=0
    }
    print chrom"\t"start"\t"end > fname
   sum+=weighted
}' intervals_with_counts.bed
   
# Rename each output file
for i in *chunk_*.bed;do
  # Pad output chunk names
  n=$(basename "$i" | sed -E 's/^chunk_([0-9]+)\.bed/\1/')
  pad=$(printf "%05d" "$n")

  # compute hash of this chunk’s contents
  hash=$(md5sum "$i" | awk '{print $1}')

  out="_${pad}${hash}.bed"
  cut -f1-4 "$i" > "$out"

  # report size
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$i")
  echo "${out}: ${bases} bases"

  rm $i
done
