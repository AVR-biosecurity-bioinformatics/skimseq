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
# $15 = max_missing


# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)


# Subset to biallelic SNPS-only
gatk SelectVariants \
	--verbosity ERROR \
	-V $2 \
	-select-type SNP \
	--restrict-alleles-to BIALLELIC \
	-O snps.vcf.gz

# Extract SNP quality scores pre filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V snps.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
	-F AF -F ExcessHet \
	-O snps.table

pigz -p${1} snps.table

# Hard-filter SNPs
if [[ ${14} == "none" ]]; then
	# use individual parameters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V snps.vcf.gz \
		-filter "QD < ${3}" --filter-name "QD${3}" \
		-filter "QUAL < ${4}" --filter-name "QUAL${4}" \
		-filter "SOR > ${5}" --filter-name "SOR${5}" \
		-filter "FS > ${6}" --filter-name "FS${6}" \
		-filter "MQ < ${7}" --filter-name "MQ${7}" \
		-filter "MQRankSum < ${8}" --filter-name "MQRankSum${8}" \
		-filter "ReadPosRankSum < ${9}" --filter-name "ReadPosRankSum${9}" \
		-filter "AF < ${10}" --filter-name "MAF${10}" \
		-filter "ExcessHet > ${11}" --filter-name "ExcessHet" \
		-filter "DP < ${12}" --filter-name "DPmin" \
		-filter "DP > ${13}" --filter-name "DPmax" \
		-O snps_tmp.vcf.gz
else
	# use custom filters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V snps.vcf.gz \
		${14} \
		-O snps_tmp.vcf.gz
fi

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V snps_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction ${15} \
	--exclude-filtered \
	-O snps_filtered_tmp.vcf.gz 

# Sort SNP vcf
gatk SortVcf \
    -I snps_filtered_tmp.vcf.gz \
    -O snps_filtered.vcf.gz 

# Extract SNP quality scores post filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V snps_filtered.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
	-F AF -F ExcessHet \
	-O snps_filtered.table

pigz -p${1} snps_filtered.table