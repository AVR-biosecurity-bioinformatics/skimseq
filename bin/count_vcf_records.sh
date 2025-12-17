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
# $8 = min_interval_gap

#GAP_BP="${8}"

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

# expand gvcf blocks to contained intervals and create bed file
# Then merge and sum number of records
# then keep only those in included intervals
#bcftools view -R included_intervals.bed -Ou ${3} \
#| bcftools convert --gvcf2vcf --fasta-ref ${4} -Ou \
#| bcftools query -f '%CHROM\t%POS0\t%POS\t1\n' \
#| {
#    IFS= read -r first || true
#    if [[ -z "${first:-}" ]]; then
#      : > ${7}.counts.bed   # no records = empty output (avoid merge error)
#    else
#      { printf '%s\n' "$first"; cat; } \
#      | bedtools merge -i - -d "$GAP_BP" -c 4 -o sum > ${7}.counts.bed
#    fi
#  }

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

# TODO: Could add callability filter here, to find sites that are covered by N reads etc