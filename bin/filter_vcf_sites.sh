#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = variant_type
# $5 = qd
# $6 = qual
# $7 = sor
# $8 = fs
# $9 = mq
# $10 = mqrs
# $11 = rprs
# $12 = maf
# $13 = mac
# $14 = eh
# $15 = dp_min
# $16 = dp_max
# $17 = max_missing
# $18 = custom_flags
# $19 = mask_bed

# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Make sure mask file is sorted and unique
bedtools sort -i ${19} | uniq > vcf_masks.bed

# Index mask bed for use in filtering
gatk IndexFeatureFile \
     -I vcf_masks.bed

# Subset to variant_type
if [[ ${4} == "snp" ]]; then 
    TYPE="-select-type SNP"
    RESTRICT="-restrict-alleles-to BIALLELIC"
elif [[ ${4} == "indel" ]]; then 
    TYPE="-select-type INDEL"
    RESTRICT="-restrict-alleles-to BIALLELIC"
elif [[ ${4} == "invariant" ]]; then 
    TYPE="-select-type NO_VARIATION"
    RESTRICT=""
else
    echo "variant_type needs to be snp, indel, or invariant"
    exit 1
fi

gatk SelectVariants \
	--verbosity ERROR \
	-V ${3} \
	$TYPE \
	$RESTRICT \
	-O variants.vcf.gz

# Set up filters
filters=()
if [[ ${18} == "none" ]]; then
	# use individual parameters
    [[ "${5}"  != NA ]] && filters+=( -filter "QD < ${5}"                --filter-name QD )
    [[ "${6}"  != NA ]] && filters+=( -filter "QUAL < ${6}"              --filter-name QUAL )
    [[ "${7}"  != NA ]] && filters+=( -filter "SOR > ${7}"               --filter-name SOR )
    [[ "${8}"  != NA ]] && filters+=( -filter "FS > ${8}"                --filter-name FS )
    [[ "${9}"  != NA ]] && filters+=( -filter "MQ < ${9}"                --filter-name MQ )
    [[ "${10}"  != NA ]] && filters+=( -filter "MQRankSum < ${10}"       --filter-name MQRankSum )
    [[ "${11}"  != NA ]] && filters+=( -filter "ReadPosRankSum < ${11}"  --filter-name ReadPosRankSum )
    [[ "${12}" != NA ]] && filters+=( -filter "MAF < ${12}"              --filter-name MAF )
    [[ "${13}" != NA ]] && filters+=( -filter "MAC < ${13}"              --filter-name MAC )
    [[ "${14}" != NA ]] && filters+=( -filter "ExcessHet > ${14}"        --filter-name ExcessHet )
    [[ "${15}" != NA ]] && filters+=( -filter "DP < ${15}"               --filter-name DPmin )
    [[ "${16}" != NA ]] && filters+=( -filter "DP > ${16}"               --filter-name DPmax )
    [[ "${17}" != NA ]] && filters+=( -filter "F_MISSING > ${17}"        --filter-name F_MISSING )

else
	# use custom filters
	filters+=(${18})
fi

# Annotate filter column for sites that fail filters
gatk VariantFiltration \
	--verbosity ERROR \
	-V variants.vcf.gz \
	"${filters[@]}" \
	--mask vcf_masks.bed --mask-name "Mask" \
	-O annot_filters.vcf.gz

# Create summary table of filters
gatk VariantsToTable \
	--verbosity ERROR \
	-V annot_filters.vcf.gz \
	-F CHROM -F POS -F TYPE -F FILTER -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum \
	-F SOR -F MAF -F MAC -F ExcessHet -F F_MISSING -F NS \
	--show-filtered \
	-O ${4}_filtered.table

pigz -p${1} ${4}_filtered.table

# Exclude filtered sites from output vcf
gatk SelectVariants \
	--verbosity ERROR \
	-V annot_filters.vcf.gz \
	--exclude-filtered \
	-O filtered_tmp.vcf.gz 

# Sort site filtered vcf
gatk SortVcf \
    -I filtered_tmp.vcf.gz \
    -O ${4}_filtered.vcf.gz 


