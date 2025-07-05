#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample name
# $4 = bam file
# $5 = ref genome
# $6 = interval hash
# $7 = interval_bed
# $8 = interval_padding
# $9 = exclude_bed
# $10 = exclude_padding

# Create list of bams to be processed
echo ${4} | tr ' ' '\n' > bam.list

# Merge all masks
touch merged_masks.bed
while read mask; do
  cat $mask | cut -f1-4 >> merged_masks.bed
done < <(echo ${9} | tr ' ' '\n')

# Merge any overlapping masks
bedtools sort -i merged_masks.bed > merged_masks_sorted.bed
bedtools merge -i merged_masks_sorted.bed -c 4 -o distinct > merged_masks_distinct.bed


# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" HaplotypeCaller \
    -R $5 \
    -I bam.list \
    -L $7 \
    -O ${3}.${6}.g.vcf.gz \
    --native-pair-hmm-threads ${1} \
    --min-base-quality-score 15 \
    --min-pruning 0 \
    --interval-padding ${8} \
    --exclude-intervals merged_masks_distinct.bed \
    --interval-exclusion-padding ${10} \
    -ERC GVCF