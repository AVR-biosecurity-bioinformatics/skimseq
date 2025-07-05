#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_bed
# $3 = excluded_bed
# $4 = Reference_genome

# Included intervals is just the reference genome bed if specific intervals were not provided
# Just keep the important columns
cat ${2} | cut -f1-3 | sed 's/\s*$/\tPASS/' > included_intervals.bed

# Merge all masks
touch merged_masks.bed
while read mask; do
  cat $mask | cut -f1-4 >> merged_exclude_masks.bed
done < <(echo ${3} | tr ' ' '\n')

# Merge any overlapping intervals
bedtools sort -i merged_exclude_masks.bed > merged_exclude_masks_sorted.bed
bedtools merge -i merged_exclude_masks_sorted.bed -c 4 -o distinct > all_intervals.bed

# Add the included intervals that arent in the mask
bedtools subtract -a included_intervals.bed -b all_intervals.bed | cut -f1-4 >> all_intervals.bed
bedtools sort -i all_intervals.bed > all_intervals_sorted.bed

# Add any remaining intervals of the genome that werent in the exclude or include masks
bedtools complement -i all_intervals_sorted.bed -g ${4}.fai | sed 's/\s*$/\tExcluded/' >> all_intervals_sorted.bed
bedtools sort -i all_intervals_sorted.bed > mask_summary.bed

# Tabulate by annotation type
awk '{len=$3-$2; sum[$4]=sum[$4]+len} END {for (anno in sum) print anno, sum[anno]}'  mask_summary.bed > mask_summary.txt
