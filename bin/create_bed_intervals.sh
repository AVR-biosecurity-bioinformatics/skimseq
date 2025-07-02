#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_splitn
# $3 = interval_n
# $4 = included_bed
# $5 = included_padding
# $6 = excluded_bed
# $7 = excluded_padding
# $8 = mitochondrial_contig
# $9 = Reference_genome

# Included intervals is just the reference genome bed if specific intervals were not provided
cat ${4} > included_intervals.bed

# Convert any scientific notation to integers
interval_splitlength=$(awk -v x="${2}" 'BEGIN {printf("%d\n",x)}') #TODO REMOVE INTERVAL SIZE
interval_n=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')

# Annotate any exluded intervals on the output bed
# Create excluded intervals list if provided, or just an empty dummy file
if [ -s ${6} ] ; then  
  bedtools slop -i ${6} -g ${9}.fai -b ${7} > excluded_intervals.bed
  # TODO: ADD EXTRA ANNOTATIONS COLUMN
else 
  touch excluded_intervals.bed
fi

# Annotate the mitochondrial contig on the output bed if present
if [ ${8} ] ; then  
  grep -i ${8} included_intervals.bed >> excluded_intervals.bed
  # TODO: ADD EXTRA ANNOTATIONS COLUMN
else 
  touch mito_contig.bed
fi

# TODO: Handle any extra filters (coverage, masks, etc)

# Detect any strings of N bases in the reference genome
java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
      REFERENCE=${9}\
      OUTPUT_TYPE=BOTH \
      OUTPUT=output.interval_list \
	    MAX_TO_MERGE=${interval_splitlength}

# Convert resulting interval list to bed format
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
	INPUT=output.interval_list \
	OUTPUT=intervals.bed

# Subtract any excluded intervals from the original list
bedtools subtract -a intervals.bed -b excluded_intervals.bed > intervals_excluded.bed

# Then add the excluded intervals to the end of it
cat excluded_intervals.bed >> intervals_excluded.bed
bedtools sort -i intervals_excluded.bed > interval_summary.bed

# Filter intervals bedfile to just good intervals (ACGT bases) for creating windows
awk '/ACGTmer/' interval_summary.bed > intervals_filtered.bed

# SPLIT INTERVALS into even groups
gatk SplitIntervals \
   -R ${9} \
   -L intervals_filtered.bed \
   --scatter-count ${interval_n} \
   -O $(pwd)
   
# Rename and convert each split output into a bed, and add padding
for i in *scattered.interval_list;do
  # Hashed output name
  HASH=$( printf '%s' "$i" | md5sum | awk '{print $1}' ) 

  # Convert resulting interval list to bed format
  java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
  	INPUT=$i \
  	OUTPUT=tmp.bed
	
  # Add padding to output
  bedtools slop -i tmp.bed -g ${9}.fai -b ${5} > ${HASH}.bed
  
  # remove intermediate files
  rm $i tmp.bed
done

# Remove temporary files
rm -f included_intervals.bed excluded_intervals.bed intervals_filtered.bed intervals_excluded.bed

