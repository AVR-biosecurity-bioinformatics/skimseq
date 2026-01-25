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
done < mask_beds.list

# Merge any overlapping intervals
bedtools sort -i concat_mask.bed \
  | bedtools merge -i - -c 4 -o distinct > merged_masks.bed
