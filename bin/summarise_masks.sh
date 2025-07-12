#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_bed
# $3 = excluded_bed
# $4 = Reference_genome

# Output will be a summary of all genomic positions and whether they were included in a mask or not

# Included intervals is just the reference genome bed if specific intervals were not provided
cat ${2} | cut -f1-3 | sed 's/\s*$/\tPASS/' > included_intervals.bed

# First add the masked intervals
cp merged_masks.bed all_intervals.bed

# Then add any included intervals that arent in the masks
bedtools sort -i included_intervals.bed -g ${4}.fai \
  | bedtools subtract -a stdin -b ${3} \
  | cut -f1-4 >> all_intervals.bed

# Then add any remaining intervals of the genome that werent in the exclude or include masks
bedtools sort -i all_intervals.bed -g ${4}.fai \
  | bedtools complement -i stdin -g ${4}.fai \
  | sed 's/\s*$/\tExcluded/' >> all_intervals.bed

# Create sorted mask summary file
bedtools sort -i all_intervals.bed -g ${4}.fai > mask_summary.bed

# Tabulate by annotation type
awk '{len=$3-$2; sum[$4]=sum[$4]+len} END {for (anno in sum) print anno, sum[anno]}'  mask_summary.bed > mask_summary.txt
