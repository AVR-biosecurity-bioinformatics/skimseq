#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = snp_qd
# $4 = snp_qual
# $5 = snp_sor
# $6 = snp_fs
# $7 = snp_mq
# $8 = snp_mqrs
# $9 = snp_rprs
# $10 = snp_maf
# $11 = snp_eh
# $12 = snp_dp_min
# $13 = snp_dp_max
# $14 = snp_custom_flags
# $15 = max_nocall
# $16 = max_missing
# $17 = gt_qual
# $18 = gt_dp_min
# $19 = gt_dp_max
# $20 = mask_bed


# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Make sure mask file is sorted and unique
bedtools sort -i ${20} | uniq > vcf_masks.bed

# Index mask bed for use in filtering
gatk IndexFeatureFile \
     -I vcf_masks.bed
     
# Subset to biallelic SNPS-only
gatk SelectVariants \
	--verbosity ERROR \
	-V $2 \
	-select-type SNP \
	--restrict-alleles-to BIALLELIC \
	-O snps.vcf.gz

# Hard-filter SNPs
if [[ ${14} == "none" ]]; then
	# use individual parameters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V snps.vcf.gz \
		-filter "QD < ${3}" --filter-name "QD" \
		-filter "QUAL < ${4}" --filter-name "QUAL" \
		-filter "SOR > ${5}" --filter-name "SOR" \
		-filter "FS > ${6}" --filter-name "FS" \
		-filter "MQ < ${7}" --filter-name "MQ" \
		-filter "MQRankSum < ${8}" --filter-name "MQRankSum" \
		-filter "ReadPosRankSum < ${9}" --filter-name "ReadPosRankSum" \
		-filter "AF < ${10}" --filter-name "AF" \
		-filter "ExcessHet > ${11}" --filter-name "ExcessHet" \
		-filter "DP < ${12}" --filter-name "DPmin" \
		-filter "DP > ${13}" --filter-name "DPmax" \
		-filter "F_MISSING > ${15}" --filter-name "F_MISSING" \
		-G-filter "GQ > ${17}" --genotype-filter-name "GQ" \
		-G-filter "DP < ${18}" --genotype-filter-name "gtDPmin" \
		-G-filter "DP > ${19}" --genotype-filter-name "gtDPmax" \
		--mask vcf_masks.bed --mask-name "Mask" \
		-O snps_tmp.vcf.gz
else
	# use custom filters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V snps.vcf.gz \
		${14} \
		--mask vcf_masks.bed --mask-name "Mask" \
		-O snps_tmp.vcf.gz
fi

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V snps_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction ${16} \
	--exclude-filtered \
	-O snps_filtered_tmp.vcf.gz 

# Sort SNP vcf
gatk SortVcf \
    -I snps_filtered_tmp.vcf.gz \
    -O snps_filtered.vcf.gz 

# Create SNP filters summary  table
gatk VariantsToTable \
	--verbosity ERROR \
	-V snps_tmp.vcf.gz \
	-F CHROM -F POS -F TYPE -F FILTER -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum \
	-F SOR -F AF -F ExcessHet -F F_MISSING -F NS \
	--show-filtered \
	-O snps_filtered.table

pigz -p${1} snps_filtered.table