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
 | bedtools sort -i stdin -g ${4} > included_intervals.bed

# Count number of VCF records overlapping intervals
#bcftools query -f '%CHROM\t%POS0\t%POS\n' ${3} \
#    | bedtools intersect -a included_intervals.bed -b - -c > ${7}.counts.bed

bedtools intersect \
    -a included_intervals.bed \
    -b ${3} \
    -sorted \
    -c > ${7}.counts.bed