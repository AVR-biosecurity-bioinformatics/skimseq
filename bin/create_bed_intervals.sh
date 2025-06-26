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
 
# Split large contigs into windows of interval_sizes
bedtools makewindows -b intervals_filtered.bed -w ${2} > intervals_filtered_windows.bed

# Split bed into multiple files of interval_n
bedtools split -i intervals_filtered_windows.bed -n ${3} -p interval

# Pad and rename output intervals to 
for i in interval.*bed;do
  echo $i
  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Add padding to output
  bedtools slop -i $i -g ${9}.fai -b ${5} > ${HASH}.bed
  
  # remove intermediate interval
  rm $i
done

# Remove temporary files
#rm included_intervals.bed excluded_intervals.bed intervals_filtered.bed

