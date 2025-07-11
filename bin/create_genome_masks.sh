#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_bed
# $3 = excluded_bed
# $4 = excluded_padding
# $5 = use_reference_hardmasks
# $6 = use_reference_softmasks
# $7 = Reference_genome

# Included intervals is just the reference genome bed if specific intervals were not provided
# Just keep the important columns
cat ${2} | cut -f1-3 > included_intervals.bed

# Create mask bed file
touch genome_masks.bed

# Add any excluded intervals to hard masked bed if file is not empty
if [ -s ${3} ] ; then  
  # Keep just the important columns
  cat ${3} \
    | cut -f1-3 \
    | sed 's/\s*$/\tExcluded/' >> genome_masks.bed
fi

# Find any N bases in reference genome
java -jar $EBROOTPICARD/picard.jar ScatterIntervalsByNs \
      --REFERENCE ${7} \
      --OUTPUT_TYPE N \
      --OUTPUT N_bases.interval_list \
  	  --MAX_TO_MERGE 1
  	    
# Convert resulting interval list to bed format 
java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    	--INPUT N_bases.interval_list \
    	--OUTPUT N_bases.bed

# If use_reference_hardmasks is true, add N base coordinates to the mask file
if [ ${5} == "true" ] ; then
  # Subset to just those inside the included intervals and add to masks bed
  bedtools subtract -a N_bases.bed -b genome_masks.bed \
    | bedtools intersect -wa -a stdin -b included_intervals.bed \
    | cut -f1-4 \
    | sed 's/Nmer/NRef/g' >> genome_masks.bed
fi

# If exclude_reference_genome_softmasks is true, add any lowercase reference genome bases to soft mask bed
if [ ${6} == "true" ] ; then
  
  # Create temporary reference genome where soft-masked regions are converted to N
  cat ${7} \
    | sed '/^[^>]/s/[^ATGC]/N/g' > tmp.fa
  
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
  	
  # Remove temporary fasta
  rm -f tmp.fa*
  	
  # Subtract any intervals that were originally N bases or already contained in masks bed
  bedtools subtract -a all_masked_bases.bed -b N_bases.bed \
    | bedtools subtract -a stdin -b genome_masks.bed \
    | bedtools intersect -wa -a stdin -b included_intervals.bed \
    | cut -f1-4 \
    | sed 's/Nmer/SoftMaskRef/g' >> genome_masks.bed
fi
