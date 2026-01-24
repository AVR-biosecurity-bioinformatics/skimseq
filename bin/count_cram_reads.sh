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

tmp_bed="${4}.counts.bed"

# count per-base depths
# Then exclude any zero counts with awk and create bed, merging abutting intervals
samtools depth \
  -@ ${1} \
  -q ${7} \
  -Q ${8} \
  ${FLAGS} \
  --reference ${5} \
  ${3} \
  |	awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' \
  | bedtools merge -i -  -c 4 -o sum > "$tmp_bed"

# bgzip output and create  tabix index
bgzip -c "$tmp_bed" > "${4}.counts.bed.gz"
tabix -f -p bed "${4}.counts.bed.gz"

# Cleanup
rm "$tmp_bed"