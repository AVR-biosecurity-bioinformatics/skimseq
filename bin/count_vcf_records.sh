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
# $8 = min_interval_gap

GAP_BP="${8}"

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

# expand gvcf blocks to contained intervals and create bed file
# Then merge and sum number of records
# then keep only those in included intervals
bcftools view -R included_intervals.bed -Ou ${3} \
| bcftools convert --gvcf2vcf --fasta-ref ${4} -Ou \
| bcftools query -f '%CHROM\t%POS0\t%POS\t1\n' \
| bedtools merge -i stdin -d "$GAP_BP" -c 4 -o sum > ${7}.counts.bed
