#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = interval_size
# $5 = counts_files

TARGET_BASES=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')
OUTDIR=$(pwd)

# Combine coverage files for all samples
bedtools unionbedg -i *counts.bed -filler 0 > combined_counts.bed

# get mean aligned bases per window
awk '{
    sum=0;
    n=0;
    for(i=4;i<=NF;i++){
        sum+=$i;
        n++
    }
    mean = (n>0 ? sum/n : 0)
    print $1"\t"$2"\t"$3"\t"mean
}' combined_counts.bed > intervals_with_depth.bed

# Use greedy algorithm to assign intervals to chunks.
# Once they reach TARGET_BASES, begin a new chunk.
awk -v target="$TARGET_BASES" -v outdir="$OUTDIR" '
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
}' intervals_with_depth.bed
   
# Rename each output to a hash
for i in *chunk_*.bed;do
  # Hashed output name
  HASH=$( md5sum "$i" | awk '{print $1}' ) 

  #adding extra character at start to ensure that other temp beds dont accidentally get passed to next process
  cut -f1-4 "$i" > _${HASH}.bed

  sum=$(awk '{sum+=$3-$2} END{print sum}' "$i")
  echo "$HASH: $sum bases"
  rm $i
done
