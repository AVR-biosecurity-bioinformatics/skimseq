#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_bed
# $3 = excluded_bed
# $4 = excluded_padding
# $5 = use_reference_hardmasks
# $6 = use_reference_softmasks
# $7 = Reference_genome

ref="${7}"

# Included intervals is just the reference genome bed if specific intervals were not provided
# Just keep the important columns
cat ${2} | cut -f1-3 > included_intervals.bed

# Create mask bed file
touch concat_masks.bed

# Add any excluded intervals to genome mask bed if file is not empty
if [ -s ${3} ] ; then  
  # Keep just the important columns
  cat ${3} \
    | cut -f1-3 \
    | sed 's/\s*$/\tExcluded/' >> concat_masks.bed
fi

# If use_reference_hardmasks is true, add N base coordinates to the mask file
if [ ${5} == "true" ] ; then

  seqkit locate -P -r -p "N" --bed \
    --id-regexp "^(\\S+)" \
    "$ref" \
  | cut -f1,2,3 \
  | bedtools merge -i - \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"NRef"}' \
  >> concat_masks.bed
fi

# If exclude_reference_genome_softmasks is true, add any lowercase reference genome bases to soft mask bed
if [ ${6} == "true" ] ; then
    seqkit locate -P -r -p "[a-z]" --bed \
    --id-regexp "^(\\S+)" \
    "$ref" \
  | cut -f1,2,3 \
  | bedtools merge -i - \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"SoftMaskRef"}' \
  >> concat_masks.bed

fi

# Restrict masks to included intervals, then sort+merge with labels
bedtools intersect -wa -a concat_masks.bed -b included_intervals.bed \
  | bedtools sort -i - \
  | bedtools merge -i - -c 4 -o distinct \
  > genome_masks.bed