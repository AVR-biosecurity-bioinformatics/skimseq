#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bin counts files
# $4 = ref genome
# $5 = binned_bed
# $6 = annotated_bins

# Create list of bams to be processed
echo ${3} | tr ' ' '\n' > count.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" FilterIntervals \
    -I count.list \
    -L ${5} \
    --annotated-intervals ${6} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output bin_filtered.interval_list \
    --extreme-count-filter-maximum-percentile 99.0 \
    --extreme-count-filter-minimum-percentile 1.0 \
    --extreme-count-filter-percentage-of-samples 90.0 \
    --low-count-filter-count-threshold 5 \
    --low-count-filter-percentage-of-samples 90.0 \
    --maximum-gc-content 0.9 \
    --minimum-gc-content 0.1 
    
# Filters to add in future:
# --maximum-mappability
# --minimum-mappability    
    
# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT bin_filtered.interval_list \
  	--OUTPUT bin_filtered.bed 
  	
# TODO: need to get the inverse of these 
  # Subset to just those inside the included intervals and add to masks bed
  bedtools intersect -wa -a tmp.bed -b included_intervals.bed | cut -f1-4 | sed 's/Nmer/SoftMaskRef/g' >> masks.bed