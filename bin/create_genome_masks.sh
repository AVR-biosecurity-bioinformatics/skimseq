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
touch mask.bed

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
  # Subset to just those inside the included intervals and add to hard masked bed
  bedtools intersect -wa -a N_bases.bed -b included_intervals.bed | cut -f1-4 | sed 's/Nmer/NRef/g' >> masks.bed
fi

# If exclude_reference_genome_softmasks is true, add any lowercase reference genome bases to soft mask bed
if [ ${6} == "true" ] ; then
  
  # Create temporary reference genome where soft-masked regions are converted to N
  cat ${2} | sed '/^[^>]/s/[^ATGC]/N/g' > tmp.fa
  
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
  	
  # Subtract any intervals that were originally hard masked bases
  bedtools subtract -a all_masked_bases.bed -b N_bases.bed > soft_masked_bases_only.bed
  
  # Subset to just those inside the included intervals and add to soft masked bed
  bedtools intersect -wa -a soft_masked_bases_only.bed -b included_intervals.bed | cut -f1-4 | sed 's/Nmer/SoftMaskRef/g' >> masks.bed

fi

# Add any excluded intervals + padding to hard masked bed
if [ -s ${3} ] ; then  
  # Keep just the important columns
  cat ${3} | cut -f1-3 > tmp.bed 
  bedtools slop -i tmp.bed -g ${6}.fai -b ${4} | sed 's/\s*$/\tExcluded/' >> masks.bed
fi


# Add mitochondrial contig to the soft masked bed
# TODO: This needs to be done in nextflow logic
#if [ ${5} ] ; then  
#  grep -i ${5} ${6}.fai |  awk '{print $1"\t0\t"$2"\tMitoContig"}' >> soft_masked.bed
#fi

# Create a summary bed listing the of proportion of genome contained within different masks:
# TODO: This is done in a new summarise_mask process
#cat soft_masked.bed > merged_masks.bed
#cat hard_masked.bed >> merged_masks.bed

#bedtools sort -i merged_masks.bed > merged_masks_sorted.bed

#bedtools complement -i merged_masks_sorted.bed -g ${6}.fai | sed 's/\s*$/\tIncluded/' >> merged_masks.bed
#bedtools sort -i merged_masks.bed > mask_summary.bed

# Tabulate by annotation type
# awk '{len=$3-$2; sum[$4]=sum[$4]+len} END {for (anno in sum) print anno, sum[anno]}' mask_summary.bed

