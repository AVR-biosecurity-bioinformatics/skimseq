#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = excluded_bed

# Concatenate any seperate mask bed files
touch concat_mask.bed
while read mask; do
  cat $mask | cut -f1-4 >> concat_mask.bed
done < <(echo ${2} | tr ' ' '\n')

# Merge any overlapping intervals
bedtools sort -i concat_mask.bed \
  | bedtools merge -i stdin -c 4 -o distinct > merged_masks.bed
