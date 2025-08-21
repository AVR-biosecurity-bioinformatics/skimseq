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

# Exclude any intervals if exclusion files are not empty
bedtools subtract -a ${5} -b ${6} > included_intervals.bed

# Count number of VCF records overlapping intervals
bedtools intersect \
    -a included_intervals.bed \
    -b ${3} \
    -c > ${7}.counts.bed
