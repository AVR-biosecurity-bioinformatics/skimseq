#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

# Find bases with a callable genotype, expanding gvcf blocks
bcftools query -R included_intervals.bed -f '%CHROM\t%POS0\t%POS\t%INFO/END\n' ${3} \
	  | awk -v OFS="\t" '
		  {
			chrom=$1; start=$2; pos=$3; endtag=$4

			# BED end: use END if present else POS (single base)
			end = (endtag=="." ? pos : endtag)
			print chrom, start, end
		  }' \
	  | bedtools merge -i - \
      > ${7}.counts.bed

# TODO: Could add callability filter here, to keep sites that are covered by N reads etc