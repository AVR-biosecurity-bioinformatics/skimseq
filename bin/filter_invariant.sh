#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = inv_dp_min
# $4 = inv_dp_max
# $5 = indel_custom_flags
# $6 = max_missing
# $7 = mask_bed

# Make sure mask file is sorted and unique
# TODO: Work out why the input mask is duplicated in the first place
bedtools sort -i ${7} | uniq > vcf_masks.bed

# Index mask bed for use in filtering
gatk IndexFeatureFile \
     -I vcf_masks.bed
     
# Subset to NO_VARIATION
gatk SelectVariants \
	--verbosity ERROR \
	-V $2 \
	-select-type NO_VARIATION \
	-O inv.vcf.gz

# Extract invariant quality scores pre filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V inv.vcf.gz \
	-F CHROM -F POS -F QUAL -F DP \
	-O inv.table

pigz -p${1} inv.table

# Hard-filter invariant
if [[ ${5} == "none" ]]; then
	# use individual parameters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V inv.vcf.gz \
		-filter "DP < ${3}" --filter-name "DPmin" \
		-filter "DP > ${4}" --filter-name "DPmax" \
		--mask vcf_masks.bed --mask-name Mask \
		-O inv_tmp.vcf.gz
else
	# use custom filters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V inv.vcf.gz \
		${5} \
		--mask vcf_masks.bed --mask-name Mask \
		-O inv_tmp.vcf.gz
fi

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V inv_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction ${6} \
	--exclude-filtered \
	-O inv_filtered_tmp.vcf.gz 

# Sort invariant vcf
gatk SortVcf \
    -I inv_filtered_tmp.vcf.gz \
    -O inv_filtered.vcf.gz 

# Extract invariant quality scores post filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V inv_filtered.vcf.gz \
	-F CHROM -F POS -F QUAL -F DP \
	-O inv_filtered.table

pigz -p${1} inv_filtered.table