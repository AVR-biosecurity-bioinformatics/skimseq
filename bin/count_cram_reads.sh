#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = cram file
# $4 = sample
# $5 = ref genome
# $6 = hc_rmdup
# $7 = hc_minbq
# $8 = hc_minmq

# Set up samtools flags
if [[ ${6} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    FLAGS="-g DUP -G UNMAP,SECONDARY,QCFAIL"
else 
    FLAGS="-G UNMAP,SECONDARY,QCFAIL,DUP"
fi

tmp_bed="tmp.bed"

# count per-base depths
samtools depth \
  -@ ${1} \
  -q ${7} \
  -Q ${8} \
  ${FLAGS} \
  --reference ${5} \
  ${3} \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}'> "$tmp_bed"
  
# Find covered tracts by merging abutting intervals and summing counts
bedtools merge -i "$tmp_bed" -c 4 -o sum \
  | bgzip -c --compress-level 9 > "${4}.covered.bed.gz"
tabix -f -p bed "${4}.covered.bed.gz"

# bgzip output and create  tabix index
bgzip -c --compress-level 9 "$tmp_bed" > "${4}.perbase.bed.gz"
tabix -f -p bed "${4}.perbase.bed.gz"

# Cleanup
rm "$tmp_bed"