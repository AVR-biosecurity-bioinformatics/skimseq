#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bam file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample
# $8 = mode
# $9 = hc_rmdup
# $10 = hc_minmq

# Set up samtools flags
if [[ ${9} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    KEEP_FLAGS="DUP"
    REMOVE_FLAGS="UNMAP,SECONDARY,QCFAIL"
else 
    REMOVE_FLAGS="UNMAP,SECONDARY,QCFAIL,DUP"
fi

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
bedtools subtract \
    -a <(cut -f1-3 "${5}") \
    -b <(cut -f1-3 "${6}") > included_intervals.bed

# Count number of aligned reads and aligned bases overlapping intervals
# Exclude by the same mapping quality and duplicate removal that will be used for gatk
samtools bedcov \
    --min-MQ ${10} \
    --reference ${4} \
    -g $KEEP_FLAGS \
    -G $REMOVE_FLAGS \
    included_intervals.bed \
    "${3}" -c > counts.bed.tmp

# Select columns based on mode
awk -v mode="${8}" 'BEGIN{OFS="\t"} 
{
    if(mode=="reads") print $1,$2,$3,$5
    else if(mode=="bases") print $1,$2,$3,$4
}' counts.bed.tmp > ${7}.counts.bed

# Remove temporary files
rm -f counts.bed.tmp
