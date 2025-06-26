#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_size
# $3 = interval_n
# $4 = included_bed
# $5 = included_padding
# $6 = excluded_bed
# $7 = excluded_padding
# $8 = mitochondrial_contig
# $9 = Reference_genome

cat ${4} > included_intervals.bed

# Convert any scientific notation to real numbers for interval_size
interval_size=$(awk -v x="${2}" 'BEGIN {printf("%d\n",x)}')
interval_n=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')

# Create excluded intervals list if provided, or just an empty dummy file
if [ -s ${6} ] ; then  
  bedtools slop -i ${6} -g ${9}.fai -b ${7} > excluded_intervals.bed
else 
  touch excluded_intervals.bed
fi

# exclude mitochondrial contig if its in the input
if [ ${8} ] ; then  
  grep -i ${8} included_intervals.bed >> excluded_intervals.bed
fi

# Subtract any excluded intervals from included list
bedtools subtract -a included_intervals.bed -b excluded_intervals.bed > intervals_filtered.bed
 
# Create interval groups

if [ -n "$interval_size" ] && [ -z "$interval_n" ]; then
    # If only interval_size is provided, split large chromosomes into intervals of size $interval_size
    # and then split into approximately $interval_size sized windows
    bedtools makewindows -b intervals_filtered.bed -w $interval_size > intervals_filtered_windows.bed

    # Calculate the total length of the genome after splitting
    total_length=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' intervals_filtered_windows.bed)

    # Calculate the number of splits based on the total length and the interval size
    n_splits=$(awk -v total_length=$total_length -v interval_size=$interval_size 'BEGIN { print ( total_length / interval_size ) }')

    # Split the BED file into the calculated number of splits
    bedtools split -i intervals_filtered_windows.bed -n $n_splits -p _split

elif [ -z "$interval_size" ] && [ -n "$interval_n" ]; then
    # If only interval_n is provided, calculate the optimal number of splits
    total_length=$(awk '{sum+=$3-$2} END {print sum}' intervals_filtered.bed)

    # Calculate the average window size based on the total length and the number of intervals
    average_window_size=$((total_length / interval_n))

    # Split the genome into windows of the calculated average window size
    bedtools makewindows -b intervals_filtered.bed -w $average_window_size > intervals_filtered_windows.bed

    # Split the BED file into the specified number of splits
    bedtools split -i intervals_filtered_windows.bed -n $interval_n -p _split

else
    # If both interval_size and interval_n are provided, first split by interval_size, then split by interval_n
    bedtools makewindows -b intervals_filtered.bed -w $interval_size > intervals_filtered_windows.bed

    # Split the resulting intervals into the specified number of splits
    bedtools split -i intervals_filtered_windows.bed -n $interval_n -p _split
fi

# Pad and rename each group of output intervals 
for i in _split.*bed;do
  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Add padding to output
  bedtools slop -i $i -g ${9}.fai -b ${5} > ${HASH}.bed
  
  # remove intermediate interval
  rm $i
done

# Remove temporary files
rm included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_filtered_windows.bed

