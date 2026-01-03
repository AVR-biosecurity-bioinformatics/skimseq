#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = cram file
# $4 = sample
# $5 = ref genome
# $6 = include_bed     
# $7 = exclude_bed
# $8 = hc_rmdup
# $9 = hc_minbq
# $10 = hc_minmq

# Set up samtools flags
if [[ ${8} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    FLAGS="-g DUP -G UNMAP,SECONDARY,QCFAIL"
else 
    FLAGS="-G UNMAP,SECONDARY,QCFAIL,DUP"
fi

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
bedtools subtract \
    -a <(cut -f1-3 "${6}") \
    -b <(cut -f1-3 "${7}") > included_intervals.bed

tmp_bed="${4}.counts.bed"

# count per-base depths
# Then exclude any zero counts with awk and create bed
# Then merge any blocks less than GAP_BP apart and sum counts
samtools depth \
  -b included_intervals.bed \
  -@ ${1} \
  -q ${9} \
  -Q ${10} \
  ${FLAGS} \
  --reference ${5} \
  ${3} \
  |	awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > "$tmp_bed"

# bgzip output and create  tabix index
bgzip -c "$tmp_bed" > "${4}.counts.bed.gz"
tabix -f -p bed "${4}.counts.bed.gz"

# Cleanup
rm "$tmp_bed"