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

# count per-base depths
samtools depth \
  -@ ${1} \
  -q ${7} \
  -Q ${8} \
  ${FLAGS} \
  --reference ${5} \
  ${3} \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' \
  | bgzip -c --compress-level 9 > "${4}.perbase.bed.gz"
  
tabix -f -p bed "${4}.perbase.bed.gz"