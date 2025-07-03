#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_bed
# $3 = excluded_bed
# $4 = excluded_padding
# $5 = mitochondrial_contig
# $6 = Reference_genome
# $7 = exclude_reference_hardmasks
# $8 = exclude_reference_softmasks

# Included intervals is just the reference genome bed if specific intervals were not provided
cat ${2} > included_intervals.bed

# Create empty output files
touch hard_masked.bed
touch soft_masked.bed

# If exclude_reference_hardmasks is true, find any existing N bases in reference genome and add to hard masked bed
if [ ${7} == "true" ] ; then
  # Find any existing N bases in reference genome and create hard masked bed
  java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
        --REFERENCE ${6} \
        --OUTPUT_TYPE N \
        --OUTPUT N_bases.interval_list \
  	    --MAX_TO_MERGE 1
  	    
  # Convert resulting interval list to bed format and add to hard masked bed
  java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    	--INPUT N_bases.interval_list \
    	--OUTPUT N_bases.bed
  cat N_bases.bed >> hard_masked.bed
fi

# Add any excluded intervals + padding to hard masked bed
if [ -s ${3} ] ; then  
  bedtools slop -i ${3} -g ${6}.fai -b ${4} | sed 's/\s*$/\tExcluded/' >> hard_masked.bed
fi

# If exclude_reference_genome_softmasks is true, add any lowercase reference genome bases to soft mask bed
if [ ${8} == "true" ] ; then
  
  # Create temporary reference genome where soft-masked regions are converted to N
  cat ${6} | sed '/^[^>]/s/[^ATGC]/N/g' > tmp.fa
  
  # Index temporary reference genome and create GATK dictionary
  samtools faidx tmp.fa
  gatk CreateSequenceDictionary -R tmp.fa

  # Find all N bases in refernce genome (soft and hard masks now)
  java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
      --REFERENCE tmp.fa \
      --OUTPUT_TYPE N \
      --OUTPUT all_masked_bases.interval_list \
	    --MAX_TO_MERGE 1
	
	# Convert resulting interval list to bed format
	java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
  	--INPUT all_masked_bases.interval_list \
  	--OUTPUT all_masked_bases.bed
  	
  # Subtrackt any intervals that were originally N's and add to soft masked bed
  bedtools subtract -a all_masked_bases.bed -b N_bases.bed | sed 's/Nmer/SoftMask/g' >> soft_masked.bed
fi

# Add mitochondrial contig to the soft masked bed
if [ ${8} ] ; then  
  grep -i ${8} included_intervals.bed | sed 's/\s*$/\tMitochondria/' >> soft_masked.bed
fi

# TODO: Add any additional masks (coverage, paralogs, mapabillity etc) to the soft masked bed here



