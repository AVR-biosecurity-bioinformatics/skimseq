#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf


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
gatk VariantFiltration \
	--verbosity ERROR \
	-V snps.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-filter "AF < 0.05" --filter-name "MAF005" \
	-filter "ExcessHet > 54.69" --filter-name "ExcessHet" \
	-filter "DP < 6" --filter-name "DPmin" \
	-filter "DP > 1500" --filter-name "DPmax" \
	-O snps_tmp.vcf.gz

# Transform filtered genotypes to nocall and keep only those with <5% missing data
gatk SelectVariants \
	--verbosity ERROR \
	-V snps_tmp.vcf.gz \
	--set-filtered-gt-to-nocall \
	--max-nocall-fraction 0.05 \
	--exclude-filtered \
	-O snps_filtered.vcf.gz 

# Extract SNP quality scores post filtering
gatk VariantsToTable \
	--verbosity ERROR \
	-V snps_filtered.vcf.gz \
	-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
	-F AF -F ExcessHet \
	-O snps_filtered.table

pigz -p${1} snps_filtered.table