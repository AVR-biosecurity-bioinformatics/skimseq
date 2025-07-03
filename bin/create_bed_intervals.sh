#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_n
# $3 = interval_size
# $4 = interval_bed     
# $5 = hard_masks_bed
# $6 = soft_masks_bed
# $7 = interval_include_hard_masks
# $8 = interval_include_soft_masks
# $9 = subdivide_intervals
# $10 = Reference_genome

# parse subdivide_intervals options
if [[ ${9} == "true" ]];   then SUBDIVISION_MODE="--subdivision-mode INTERVAL_SUBDIVISION "; \
else SUBDIVISION_MODE="--subdivision-mode  BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"; fi

# Convert any scientific notation to integers
interval_size=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')
interval_n=$(awk -v x="${2}" 'BEGIN {printf("%d\n",x)}')

# included_intervals is just reference genome bed if specific intervals were not provided
cat ${4} > included_intervals.bed

# Apply hard masks if interval_include_hard_masks is false (default)
if [ ${7} == "false" ] ; then
   # Subtract any of the excluded intervals
  bedtools subtract -a ${4} -b ${5} > intervals_filtered.bed
  mv intervals_filtered.bed included_intervals.bed
fi

# Apply soft masks if interval_include_soft_masks is false (default)
if [ ${8} == "false" ] ; then
   # Subtract any of the excluded intervals
  bedtools subtract -a ${4} -b ${6} > intervals_filtered.bed
  mv intervals_filtered.bed included_intervals.bed
fi

# Calculate number of groups
if [ "$interval_size" -ge 0 ] && [ "$interval_n" -eq -1 ]; then
    # Only interval_size provided
    
    # Calculate the total length and number of splits based on interval_size
    total_length=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' included_intervals.bed)
    n_splits=$(awk -v total_length=$total_length -v interval_size=$interval_size 'BEGIN { print ( total_length / interval_size ) }')
    
elif [ "$interval_size" -eq -1 ] && [ "$interval_n" -ge 0 ]; then
    # Only interval_n provided
    n_splits=$interval_n
    
elif [ "$interval_size"  -ge 0 ] && [ "$interval_n" -ge 0 ]; then
    # Both interval_size and interval_n provided, prefer interval_n
    n_splits=$interval_n
    
else
    # If neither are provided, dont do any splitting
    n_splits=1
fi

# Group current intervals into even groups of approximately even base content
# Optionally subdivide furthr
gatk SplitIntervals \
   -R ${10} \
   -L included_intervals.bed \
   --scatter-count ${n_splits} \
   -O $(pwd) \
   $SUBDIVISION_MODE
   
# Rename and convert each split output into a bed, and add padding
for i in *scattered.interval_list;do
  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Convert resulting interval list to bed format
  java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
  	--INPUT $i \
  	--OUTPUT tmp.bed
	
	#adding extra character at start to ensure that other temp beds dont accidentally get passed to next process
	cat tmp.bed | cut -f1-4 > _${HASH}.bed
	
  # remove intermediate files
  rm $i tmp.bed
done

# Remove temporary bed files
rm -f intervals_filtered.bed included_intervals.bed

