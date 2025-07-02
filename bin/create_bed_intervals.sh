#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_n
# $3 = interval_break_n_length
# $4 = subdivide_intervals
# $5 = included_bed
# $6 = excluded_bed
# $7 = excluded_padding
# $8 = mitochondrial_contig
# $9 = Reference_genome
# $9 = interval_break_n

# parse subdivide_intervals options
if [[ ${4} == "true" ]];   then SUBDIVISION_MODE="--subdivision-mode INTERVAL_SUBDIVISION "; \
else SUBDIVISION_MODE="--subdivision-mode  BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"; fi

# Included intervals is just the reference genome bed if specific intervals were not provided
cat ${5} > included_intervals.bed

# Annotate any exluded intervals on the output bed
# Create excluded intervals list if provided, or just an empty dummy file
if [ -s ${6} ] ; then  
  bedtools slop -i ${6} -g ${9}.fai -b ${7} | sed 's/\s*$/\tExcluded/' > excluded_intervals.bed
else 
  touch excluded_intervals.bed
fi

# Add mitochondrial contig to the exclusion list
if [ ${8} ] ; then  
  grep -i ${8} included_intervals.bed | sed 's/\s*$/\tMitochondria/' >> excluded_intervals.bed
fi

# TODO: Add any extra filters (coverage, masks, etc) to the exclusion list here

# Locate appropriate breakpoints in the genome using strings of Ns
if [ ${4} == "true" ] ; then
  # Detect any strings of N bases in the reference genome to define breakpoints
  java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
      --REFERENCE ${9}\
      --OUTPUT_TYPE BOTH \
      --OUTPUT breakpoints.interval_list \
	    --MAX_TO_MERGE ${3}
	    
  # Convert resulting interval list to bed format
  java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
  	--INPUT breakpoints.interval_list \
  	--OUTPUT tmp.bed
  	
  # Subset the breakpoints bed to just the included intervals (this is just the genome bed if not provided)
  bedtools intersect -wa -a tmp.bed -b included_intervals.bed | cut -f1-4 > breakpoints.bed
  
  # Subtract any of the excluded intervals - and subset  to just genotypable intervals (ACGTmers ) 
  bedtools subtract -a breakpoints.bed -b excluded_intervals.bed | awk '/ACGTmer/' > intervals_filtered.bed
else
  # Subtract any of the excluded intervals - and make summary file
  bedtools subtract -a included_intervals.bed -b excluded_intervals.bed > intervals_filtered.bed
fi


# SPLIT INTERVALS into even groups
gatk SplitIntervals \
   -R ${9} \
   -L intervals_filtered.bed \
   --scatter-count ${2} \
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
	
	cat tmp.bed | cut -f1-4 > interval_${HASH}.bed
	
  # remove intermediate files
  rm $i tmp.bed
done

# Remove temporary bed files
rm -f included_intervals.bed excluded_intervals.bed intervals_filtered.bed breakpoints.bed breakpoints_excluded.bed breakpoints.interval_list

