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

# Add any included intervals that arent in the mask to the summary file
bedtools subtract -a included_intervals.bed -b cat ${3} \
  | bedtools sort -i stdin \
  | cut -f1-4 > all_intervals.bed

# Add any remaining intervals of the genome that werent in the exclude or include masks to summary file
bedtools complement -i all_intervals_sorted.bed -g ${4}.fai \
  | bedtools sort -i stdin \
  | sed 's/\s*$/\tExcluded/' >> all_intervals.bed

bedtools sort -i all_intervals.bed > mask_summary.bed

# Tabulate by annotation type
awk '{len=$3-$2; sum[$4]=sum[$4]+len} END {for (anno in sum) print anno, sum[anno]}'  mask_summary.bed > mask_summary.txt
