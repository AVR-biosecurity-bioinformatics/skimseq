#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = reference_genome
# $3 = interval_size
# $4 = included_intervals
# $5 = included_padding
# $6 = excluded_intervals
# $7 = excluded_padding
# $8 = mitochondrial contig

# Use whole genome as the included intervals, unless interval list is provided
if [ ${4} ] ; then  
  cat ${4} > included_intervals.bed
else 
  awk '{print $1"\t0\t"$2}' ${2} > included_intervals.bed
fi

# Create excluded intervals list if provided
if [ ${6} ] ; then  
  bedtools slop -i ${6} -b ${7} > excluded_intervals.bed
else 
  touch excluded_intervals.bed
fi

# exclude mitochondrial contig if its in the input
if [ ${8} ] ; then  
  grep -i -f ${4} -i included_intervals.bed >> excluded_intervals.bed
fi

# Subtract any excluded intervals from included list
bedtools subtract -a included_intervals.bed -b excluded_intervals.bed > intervals_filtered.bed
 
# Split large contigs into windows of interval_sizes
bedtools makewindows -b intervals_filtered.bed -w ${3} > intervals_filtered_windows.bed

# Calculate number of chunks
total_length=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  intervals_filtered_windows.bed)
n_splits=$(awk -v total_length=$total_length -v interval_size=${3} 'BEGIN { print  ( total_length / interval_size ) }')

# Split bed into multiple files of approximately interval size
bedtools split -i intervals_filtered_windows.bed -n $n_splits -p interval

# Pad and rename output intervals to 
for i in *interval.*bed

  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Add padding to output
  bedtools slop -i $i -b ${5} > ${HASH}.bed
  
  # remove intermediate interval
  rm $i
fi

# Remove temporary files
rm included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_filtered_windows.bed

