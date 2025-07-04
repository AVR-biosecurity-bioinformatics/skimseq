#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = bin_size
# $3 = include_bed     
# $4 = Reference_genome

gatk PreprocessIntervals \
    -R ${3} \
    -L ${3} \
    --bin-length ${2} \
    --padding 0 \
    -O binned_intervals.interval_list

# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT binned_intervals.interval_list \
  	--OUTPUT tmp.bed
  	
cat tmp.bed | cut -f1-4 > binned_intervals.bed