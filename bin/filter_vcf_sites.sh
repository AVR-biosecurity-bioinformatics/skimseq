#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = variant_type
# $5 = mask_bed
# $6 = flags

# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Capture fixed args
CPUS=$1
MEM=$2
VCF=$3
VTYPE=$4
MASK=$5

# Capture everything after input 5 as a filter arg
shift 5
FILTER_ARGS=( "$@" )

echo "FILTER_ARGS count: ${#FILTER_ARGS[@]}"
for i in "${!FILTER_ARGS[@]}"; do printf 'FILTER_ARGS[%d]=%q\n' "$i" "${FILTER_ARGS[$i]}"; done

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${MEM} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# Make sure mask file is sorted and unique
cat $MASK > vcf_masks.bed

# Index mask bed for use in filtering
gatk IndexFeatureFile \
     -I vcf_masks.bed

# Subset to variant_type
TYPE= RESTRICT=
case "${VTYPE}" in
  snp)       TYPE="-select-type SNP";        RESTRICT="-restrict-alleles-to BIALLELIC" ;;
  indel)     TYPE="-select-type INDEL";      RESTRICT="-restrict-alleles-to BIALLELIC" ;;
  invariant) TYPE="-select-type NO_VARIATION" ;;
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" SelectVariants \
	--verbosity ERROR \
	-V ${VCF} \
	$TYPE \
	$RESTRICT \
	-O tmp.vcf.gz

# Annotate filter column for sites that fail filters
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" VariantFiltration \
	--verbosity ERROR \
	-V tmp.vcf.gz \
	"${FILTER_ARGS[@]}" \
	--mask vcf_masks.bed --mask-name "Mask" \
	-O tmp_annot.vcf.gz

# Create summary table of filters
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" VariantsToTable \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	-F CHROM -F POS -F TYPE -F FILTER -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum \
	-F SOR -F MAF -F MAC -F ExcessHet -F F_MISSING -F NS \
	--show-filtered \
	-O ${VTYPE}_filtered.table

pigz -p${CPUS} ${VTYPE}_filtered.table

# Exclude filtered sites from output vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" SelectVariants \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	--exclude-filtered \
	-O tmp_filtered.vcf.gz 

# Sort site filtered vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" SortVcf \
    -I tmp_filtered.vcf.gz \
    -O ${VTYPE}_filtered.vcf.gz 

# Remove temporary vcf files
rm -f tmp*