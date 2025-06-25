#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_n
# $3 = included_bed
# $4 = included_padding
# $5 = excluded_bed
# $6 = excluded_padding
# $7 = mitochondrial_contig
# $8 = Reference_genome

cat ${3} > included_intervals.bed

# Create excluded intervals list if provided, or just an empty dummy file
if [ -s ${5} ] ; then  
  bedtools slop -i ${5} -g ${8}.fai -b ${6} > excluded_intervals.bed
else 
  touch excluded_intervals.bed
fi

# exclude mitochondrial contig if its in the input
if [ ${7} ] ; then  
  grep -i ${7} included_intervals.bed >> excluded_intervals.bed
fi

# Subtract any excluded intervals from included list
bedtools subtract -a included_intervals.bed -b excluded_intervals.bed > intervals_filtered.bed
 
 
# OLD METHOD BASED ON INTERVAL SIZE
# Split large contigs into windows of interval_sizes
#bedtools makewindows -b intervals_filtered.bed -w ${2} > intervals_filtered_windows.bed

# Calculate number of chunks
#total_length=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  intervals_filtered_windows.bed)
#n_splits=$(awk -v total_length=$total_length -v interval_size=${2} 'BEGIN { print  ( total_length / interval_size ) }')

# Split bed into multiple files of approximately interval size
bedtools split -i intervals_filtered_windows.bed -n ${2} -p interval

# Pad and rename output intervals to 
for i in interval.*bed;do
  echo $i
  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Add padding to output
  bedtools slop -i $i -g ${8}.fai -b ${4} > ${HASH}.bed
  
  # remove intermediate interval
  rm $i
done

# Remove temporary files
rm included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_filtered_windows.bed

