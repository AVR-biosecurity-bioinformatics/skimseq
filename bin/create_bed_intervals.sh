#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_n
# $3 = interval_nbreaks
# $4 = subdivide_intervals
# $5 = included_bed
# $6 = excluded_bed
# $7 = excluded_padding
# $8 = mitochondrial_contig
# $9 = Reference_genome

# parse filtering options as flags
if [[ ${4} == "true" ]];   then SUBDIVISION_MODE="--subdivision-mode INTERVAL_SUBDIVISION "; \
else SUBDIVISION_MODE="--subdivision-mode  BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"; fi

# Included intervals is just the reference genome bed if specific intervals were not provided
cat ${5} > included_intervals.bed

# Convert any scientific notation to integers
interval_n=$(awk -v x="${2}" 'BEGIN {printf("%d\n",x)}')
interval_nbreaks=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')

# Annotate any exluded intervals on the output bed
# Create excluded intervals list if provided, or just an empty dummy file
if [ -s ${6} ] ; then  
  bedtools slop -i ${6} -g ${9}.fai -b ${7} | sed 's/\s*$/\tExcluded/' > excluded_intervals.bed
else 
  touch excluded_intervals.bed
fi

# Add mitochondrial contig to the exclusion list
if [ ${8} ] ; then  
  grep -i ${8} included_intervals.bed| sed 's/\s*$/\tMitochondria/' >> excluded_intervals.bed
fi

# TODO: Add any extra filters (coverage, masks, etc) to the exclusion list here

# Detect any strings of N bases in the reference genome to define breakpoints
java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
      --REFERENCE ${9}\
      --OUTPUT_TYPE BOTH \
      --OUTPUT breakpoints.interval_list \
	    --MAX_TO_MERGE ${interval_nbreaks}

# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
	--INPUT breakpoints.interval_list \
	--OUTPUT tmp.bed
	
# Remove unnecessary columns from bed (score, direction)
cat tmp.bed | cut -f1-4 > breakpoints.bed

# Subtract any excluded intervals from the original list
bedtools subtract -a breakpoints.bed -b excluded_intervals.bed > breakpoints_excluded.bed

# Then add the excluded intervals to the end of it to create a summary
cat excluded_intervals.bed >> breakpoints_excluded.bed
bedtools sort -i breakpoints_excluded.bed > interval_summary.bed

# Filter intervals bedfile to just good intervals (ACGT bases) for creating windows
awk '/ACGTmer/' interval_summary.bed > intervals_filtered.bed

# SPLIT INTERVALS into even groups
gatk SplitIntervals \
   -R ${9} \
   -L intervals_filtered.bed \
   --scatter-count ${interval_n} \
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
	
	cat tmp.bed | cut -f1-4 > ${HASH}.bed
	
  # Add padding to output - DO THIS IN GATK
  #bedtools slop -i tmp.bed -g ${9}.fai -b ${5} > ${HASH}.bed
  
  # remove intermediate files
  rm $i tmp.bed
done

# Remove temporary files
#rm -f included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_excluded.bed

