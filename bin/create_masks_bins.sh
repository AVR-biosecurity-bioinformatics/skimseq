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
  	--OUTPUT tmp.bed 
  	
cat tmp.bed | cut -f1-3 | sed 's/\s*$/\tGCFilt/' > GC_filtered.bed

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
  	--OUTPUT tmp.bed 
  	
cat tmp.bed | cut -f1-3 | sed 's/\s*$/\tCovFilt/' > cov_filtered.bed
    
# Combine the two filters
cat GC_filtered.bed > tmp2.bed
cat cov_filtered.bed >> tmp2.bed
bedtools sort -i tmp2.bed > bin_filtered.bed

# Get those that were filtered out
bedtools subtract -a ${5} -b bin_filtered.bed | cut -f1-4 > bin_masked.bed

# Filters to add in future:
# --maximum-mappability
# --minimum-mappability    
    