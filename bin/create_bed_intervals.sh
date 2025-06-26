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

# Convert any scientific notation to integers
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
if [ "$interval_size" -ge 0 ] && [ "$interval_n" -eq -1 ]; then
    # Only interval_size provided
    # and then split into approximately $interval_size sized windows
    bedtools makewindows -b intervals_filtered.bed -w $interval_size > intervals_filtered_windows.bed

    # Calculate the total length and number of splits based on interval_size
    total_length=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' intervals_filtered_windows.bed)
    n_splits=$(awk -v total_length=$total_length -v interval_size=$interval_size 'BEGIN { print ( total_length / interval_size ) }')

    # Split the BED file into the calculated number of splits
    bedtools split -i intervals_filtered_windows.bed -n $n_splits -p _split

elif [ "$interval_size" -eq -1 ] && [ "$interval_n" -ge 0 ]; then
    # Only interval_n provided
    total_length=$(awk '{sum+=$3-$2} END {print sum}' intervals_filtered.bed)

    # Calculate the average window size and split accordingly
    average_window_size=$((total_length / interval_n))
    bedtools makewindows -b intervals_filtered.bed -w $average_window_size > intervals_filtered_windows.bed
    bedtools split -i intervals_filtered_windows.bed -n $interval_n -p _split

elif [ "$interval_size"  -ge 0 ] && [ "$interval_n" -ge 0 ]; then
    # Both interval_size and interval_n provided
    bedtools makewindows -b intervals_filtered.bed -w $interval_size > intervals_filtered_windows.bed
    bedtools split -i intervals_filtered_windows.bed -n $interval_n -p _split
else
    # If neither are provided, use the whole inerval
    cat intervals_filtered.bed > _split.0001.bed
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
rm -f included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_filtered_windows.bed

