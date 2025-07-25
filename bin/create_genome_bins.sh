#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = bin_size
# $3 = include_bed     
# $4 = Reference_genome

gatk PreprocessIntervals \
    -R ${4} \
    -L ${3} \
    --bin-length ${2} \
    --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O binned_intervals.interval_list

# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT binned_intervals.interval_list \
  	--OUTPUT tmp.bed
  	
cat tmp.bed | cut -f1-3 > binned_intervals.bed

# Annotate intervals
gatk AnnotateIntervals \
    -R ${4} \
    -L binned_intervals.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O annotated_intervals.tsv