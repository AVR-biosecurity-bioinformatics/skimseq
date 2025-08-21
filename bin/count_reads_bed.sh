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

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") > included_intervals.bed

# Count number of reads overlapping intervals
bedtools coverage \
    -a included_intervals.bed \
    -b ${3} \
    -g ${4}.fai \
    -counts > ${7}.counts.bed
