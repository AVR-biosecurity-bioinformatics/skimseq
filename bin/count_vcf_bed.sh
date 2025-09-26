#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

bedtools intersect \
    -a included_intervals.bed \
    -b ${3} \
    -sorted \
    -g ${4}.fai \
    -c > ${7}.counts.bed