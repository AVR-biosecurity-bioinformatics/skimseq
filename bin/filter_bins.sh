#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bin counts files
# $4 = ref genome
# $6 = binned_bed
# $7 = annotated_bins

# Create list of bams to be processed
echo ${3} | tr ' ' '\n' > count.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" FilterIntervals \
    -R ${4} \
    -I count.list \
    -L ${6} \
    --annotated-intervals ${7} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --format TSV \
    --output bin_filtered.interval_list 
    
# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT bin_filtered.interval_list \
  	--OUTPUT bin_filtered.bed 