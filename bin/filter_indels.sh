#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = indel_qd
# $4 = indel_qual
# $5 = indel_fs
# $6 = indel_rprs
# $7 = indel_maf
# $8 = indel_mac
# $9 = indel_eh
# $10 = indel_dp_min
# $11 = indel_dp_max
# $12 = indel_custom_flags
# $13 = max_nocall
# $14 = max_missing
# $15 = gt_qual
# $16 = gt_dp_min
# $17 = gt_dp_max
# $18 = mask_bed


# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Make sure mask file is sorted and unique
# TODO: Work out why the input mask is duplicated in the first place
bedtools sort -i ${18} | uniq > vcf_masks.bed

# Index mask bed for use in filtering
gatk IndexFeatureFile \
     -I vcf_masks.bed
 
# Subset to biallelic INDELS 
#NOTE: MIXED events (INDEL + SNP) will be lost
gatk SelectVariants \
	--verbosity ERROR \
	-V $2 \
	-select-type INDEL \
	--restrict-alleles-to BIALLELIC \
	-O indels.vcf.gz 

# Hard-filter INDELS
if [[ ${12} == "none" ]]; then
	# use individual parameters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V indels.vcf.gz \
		-filter "QD < ${3}" --filter-name "QD" \
		-filter "QUAL < ${4}" --filter-name "QUAL" \
		-filter "FS > ${5}" --filter-name "FS" \
		-filter "ReadPosRankSum < ${6}" --filter-name "ReadPosRankSum" \
		-filter "MAF < ${7}" --filter-name "MAF" \
		-filter "MAC < ${8}" --filter-name "MAC" \
		-filter "ExcessHet > ${9}" --filter-name "ExcessHet" \
		-filter "DP < ${10}" --filter-name "DPmin" \
		-filter "DP > ${11}" --filter-name "DPmax" \
		-filter "F_MISSING > ${13}" --filter-name "F_MISSING" \
		-G-filter "GQ < ${15}" --genotype-filter-name "GQ" \
		-G-filter "DP < ${16}" --genotype-filter-name "gtDPmin" \
		-G-filter "DP > ${17}" --genotype-filter-name "gtDPmax" \
		--mask vcf_masks.bed --mask-name "Mask" \
		-O indels_tmp.vcf.gz
else
	# use custom filters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V indels.vcf.gz \
		${12} \
		--mask vcf_masks.bed --mask-name "Mask" \
		-O indels_tmp.vcf.gz
fi

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V indels_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction ${14} \
	--exclude-filtered \
	-O indels_filtered_tmp.vcf.gz

# Sort INDEL vcf
gatk SortVcf \
    -I indels_filtered_tmp.vcf.gz \
    -O indels_filtered.vcf.gz 

# Create INDEL filters summary  table
gatk VariantsToTable \
	--verbosity ERROR \
	-V indels_tmp.vcf.gz \
	-F CHROM -F POS -F TYPE -F FILTER -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum \
	-F SOR -F MAF -F MAC -F ExcessHet -F F_MISSING -F NS \
	--show-filtered \
	-O indels_filtered.table

pigz -p${1} indels_filtered.table
