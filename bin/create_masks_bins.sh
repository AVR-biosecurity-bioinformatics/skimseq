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
# $7 = bin_gc_lower
# $8 = bin_gc_upper
# $9 = bin_min_reads
# $10 = bin_lower_read_perc
# $11 = bin_upper_read_perc
# $12 = bin_filter_perc_samples


# Create list of bams to be processed
echo ${3} | tr ' ' '\n' > count.list

# First apply GC percentage filters to bins
gatk --java-options "-Xmx${2}G" FilterIntervals \
    -L ${5} \
    --annotated-intervals ${6} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output gc_filtered.interval_list \
    --minimum-gc-content ${7} \
    --maximum-gc-content ${8}
    
# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT gc_filtered.interval_list \
  	--OUTPUT gc_passing.bed 
  	
# Get the ones that failed the GC filter
bedtools subtract -a ${5} -b gc_passing.bed \
  | cut -f1-3 \
  | sed 's/\s*$/\tGCFilt/' > gc_failed.bed

# Then apply coverage filters to bins
gatk --java-options "-Xmx${2}G" FilterIntervals \
    -I count.list \
    -L ${5} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output cov_filtered.interval_list \
    --low-count-filter-count-threshold ${9} \
    --low-count-filter-percentage-of-samples ${12} \
    --extreme-count-filter-minimum-percentile ${10} \
    --extreme-count-filter-maximum-percentile ${11} \
    --extreme-count-filter-percentage-of-samples ${12}


# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    --INPUT cov_filtered.interval_list \
  	--OUTPUT cov_passing.bed 
  	
# Get the ones that failed the coverage filter
bedtools subtract -a ${5} -b cov_passing.bed \
  | cut -f1-3 \
  | sed 's/\s*$/\tCovFilt/' > cov_failed.bed
    
# Output masks for those which failed both filters
cat gc_failed.bed > failed.bed
cat cov_failed.bed >> failed.bed
bedtools sort -i failed.bed > bin_masked.bed

# Output filter file for those which passed both filters
bedtools subtract -a ${5} -b bin_masked.bed | cut -f1-4 > bin_filtered.bed

# Filters to add in future:
# --maximum-mappability
# --minimum-mappability    
    