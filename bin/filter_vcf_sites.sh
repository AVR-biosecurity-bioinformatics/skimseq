#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = variant_type
# $5 = flags
# $6 = mask_bed

# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# Make sure mask file is sorted and unique
bedtools sort -i ${6} | uniq > vcf_masks.bed

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

gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SelectVariants \
	--verbosity ERROR \
	-V ${3} \
	$TYPE \
	$RESTRICT \
	-O tmp.vcf.gz

# Annotate filter column for sites that fail filters
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" VariantFiltration \
	--verbosity ERROR \
	-V tmp.vcf.gz \
	"${5}" \
	--mask vcf_masks.bed --mask-name "Mask" \
	-O tmp_annot.vcf.gz

# Create summary table of filters
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" VariantsToTable \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	-F CHROM -F POS -F TYPE -F FILTER -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum \
	-F SOR -F MAF -F MAC -F ExcessHet -F F_MISSING -F NS \
	--show-filtered \
	-O ${4}_filtered.table

pigz -p${1} ${4}_filtered.table

# Exclude filtered sites from output vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SelectVariants \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	--exclude-filtered \
	-O tmp_filtered.vcf.gz 

# Sort site filtered vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SortVcf \
    -I tmp_filtered.vcf.gz \
    -O ${4}_filtered.vcf.gz 

# Remove temporary vcf files
rm -f tmp*