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
# Just keep the important columns
cat ${2} | cut -f1-3 > included_intervals.bed

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
  	    
  # Convert resulting interval list to bed format 
  java -jar $EBROOTPICARD/picard.jar IntervalListToBed \
    	--INPUT N_bases.interval_list \
    	--OUTPUT N_bases.bed
    	
  # Subset to just those inside the included intervals and add to hard masked bed
  bedtools intersect -wa -a N_bases.bed -b included_intervals.bed | cut -f1-4 | sed 's/Nmer/NRef/g' >> hard_masked.bed
fi

# Add any excluded intervals + padding to hard masked bed
if [ -s ${3} ] ; then  
  # Keep just the important columns
  cat ${3} | cut -f1-3 > tmp.bed 
  bedtools slop -i tmp.bed -g ${6}.fai -b ${4} | sed 's/\s*$/\tExcluded/' >> hard_masked.bed
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
  	
  # Subtract any intervals that are already in the hard masked bases
  bedtools subtract -a all_masked_bases.bed -b hard_masked.bed > soft_masked_bases_only.bed
  
  # Subset to just those inside the included intervals and add to soft masked bed
  bedtools intersect -wa -a soft_masked_bases_only.bed -b included_intervals.bed | cut -f1-4 | sed 's/Nmer/SoftMaskRef/g' >> soft_masked.bed
fi

# Add mitochondrial contig to the soft masked bed
if [ ${5} ] ; then  
  grep -i ${5} ${6}.fai |  awk '{print $1"\t0\t"$2"\tMitoContig"}' >> soft_masked.bed
fi

# TODO: Add any additional masks (coverage, paralogs, mapabillity etc) to the soft masked bed here
#
#

# Create a summary bed listing the of proportion of genome contained within different masks:
cat soft_masked.bed > merged_masks.bed
cat hard_masked.bed >> merged_masks.bed

bedtools sort -i merged_masks.bed > merged_masks_sorted.bed

bedtools complement -i merged_masks_sorted.bed -g ${6}.fai | sed 's/\s*$/\tIncluded/' >> merged_masks.bed
bedtools sort -i merged_masks.bed > mask_summary.bed

