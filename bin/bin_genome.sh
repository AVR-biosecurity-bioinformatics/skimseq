#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = bin_size
# $3 = include_bed     
# $4 = Reference_genome

# Convert any scientific notation to integers
interval_n=$(awk -v x="${2}" 'BEGIN {printf("%d\n",x)}')
interval_size=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')

# included_intervals is just reference genome bed if specific intervals were not provided
cat ${3} > included_intervals.bed

gatk PreprocessIntervals \
    -R reference.fa \
    -L ${3} \
    --bin-length ${2} \
    --padding 0 \
    -O binned_intervals.interval_list

# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT binned_intervals.interval_list \
  	--OUTPUT  binned_intervals.bed