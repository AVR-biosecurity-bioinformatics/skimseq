#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = interval_size
# $4 = include_bed     
# $5 = exclude_bed
# $6 = Reference_genome
# $7 = bam_files

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

TARGET_BASES=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')
OUTDIR=$(pwd)


# Combine coverage files for all samples
bedtools unionbedg -i *counts.bed -filler 0 > combined_counts.bed

# get total aligned bases per window
awk '{
    sum=0;
    for(i=4;i<=NF;i++) sum+=$i;
    print $1"\t"$2"\t"$3"\t"sum
}' combined_counts.bed > intervals_with_depth.bed

# Use greedy algorithm to chunk
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
