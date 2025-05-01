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
# $8 = indel_eh
# $9 = indel_dp_min
# $10 = indel_dp_max
# $11 = indel_custom_flags
# $12 = max_missing


# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)


# Subset to biallelic INDELS 
#NOTE: MIXED events (INDEL + SNP) will be lost
gatk SelectVariants \
	--verbosity ERROR \
	-V $2 \
	-select-type INDEL \
	--restrict-alleles-to BIALLELIC \
	-O indels.vcf.gz 

# Extract INDEL quality scores pre filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V indels.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
	-F AF -F ExcessHet \
	-O indels.table

pigz -p${1} indels.table

# Hard-filter INDELS
if [[ ${11} == "none" ]]; then
	# use individual parameters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V indels.vcf.gz \
		-filter "QD < ${3}" --filter-name "QD${3}" \
		-filter "QUAL < ${4}" --filter-name "QUAL${4}" \
		-filter "FS > ${5}" --filter-name "FS${5}" \
		-filter "ReadPosRankSum < ${6}" --filter-name "ReadPosRankSum${6}" \
		-filter "AF < ${7}" --filter-name "MAF${7}" \
		-filter "ExcessHet > ${8}" --filter-name "ExcessHet" \
		-filter "DP < ${9}" --filter-name "DPmin" \
		-filter "DP > ${10}" --filter-name "DPmax" \
		-O indels_tmp.vcf.gz
else
	# use custom filters
	gatk VariantFiltration \
		--verbosity ERROR \
		-V indels.vcf.gz \
		${11} \
		-O indels_tmp.vcf.gz
fi

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V indels_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction ${12} \
	--exclude-filtered \
	-O indels_filtered_tmp.vcf.gz

# Sort INDEL vcf
gatk SortVcf \
    -I indels_filtered_tmp.vcf.gz \
    -O indels_filtered.vcf.gz 

# Extract INDEL quality scores post filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V indels_filtered.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
	-F AF -F ExcessHet \
	-O indels_filtered.table

pigz -p${1} indels_filtered.table 