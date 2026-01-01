#!/bin/bash
set -e
set -u

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
SAMPLE="${3}"
CRAM="${4}"
REF="${5}"
IHASH="${6}"
INTERVAL_BED="${7}"
EXCLUDE_BED="${8}"

# 1GB of memory should be retained outside the java heap
java_mem=$(( MEM_GB - 1 ))
if (( java_mem < 1 )); then
  java_mem=1
fi

# parse filtering options as flags
if [[ "${RMDUP}" == "false" ]]; then
  RMDUP="-DF NotDuplicateReadFilter"
else
  RMDUP=""
fi

# PCR-free true = disable PCR indel error model
if [[ "${PCR_FREE}" == "true" ]]; then
  PCR_FREE="--pcr-indel-model NONE"
else
  PCR_FREE=""
fi

# GATK arg is the inverse: dont-use-soft-clipped-bases
if [[ "${USE_SOFTCLIPPED_BASES}" == "true" ]]; then
  DONT_USE_SOFTCLIPPED="false"
else
  DONT_USE_SOFTCLIPPED="true"
fi

# call variants by sample * interval chunk
# NOTE: need to use assembly region padding rather than interval_padding to avoid overlapping variants
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=${CPUS}" HaplotypeCaller \
    -R "${REF}" \
    -I ${CRAM} \
    -L "${INTERVAL_BED}" \
    --native-pair-hmm-threads "${CPUS}" \
    --assembly-region-padding "${INTERVAL_PAD}" \
    --exclude-intervals "${EXCLUDE_BED}" \
    --interval-exclusion-padding "${EXCLUDE_PAD}" \
    --interval-merging-rule ALL \
    --min-pruning "${MIN_PRUNING}" \
    --min-dangling-branch-length "${MIN_DANGLE}" \
    --max-reads-per-alignment-start "${MAX_READS_STARTPOS}" \
    ${RMDUP} \
    ${PCR_FREE} \
    --min-base-quality-score "${MINBQ}" \
    --minimum-mapping-quality "${MINMQ}" \
    --read-filter AmbiguousBaseReadFilter \
    --ambig-filter-bases "${MAX_AMBIG_BASES}" \
    --read-filter FragmentLengthReadFilter \
    --min-fragment-length "${MIN_FRAGMENT_LENGTH}" \
    --max-fragment-length "${MAX_FRAGMENT_LENGTH}" \
    --dont-use-soft-clipped-bases "${DONT_USE_SOFTCLIPPED}" \
    --read-filter OverclippedReadFilter \
    --filter-too-short "${MIN_ALIGNED_LENGTH}" \
    --mapping-quality-threshold-for-genotyping "${MINMQ}" \
    --assembly-region-out "${SAMPLE}.${IHASH}.assembly.tsv" \
    -ploidy "${PLOIDY}" \
    --heterozygosity "${HET}" \
    --heterozygosity-stdev "${HET_SD}" \
    --indel-heterozygosity "${INDEL_HET}" \
    -ERC GVCF \
    -O tmp.g.vcf.gz \
    2> >(tee -a "${IHASH}.${SAMPLE}.stderr.log" >&2)


# Extract readgroups from cram for embedding in VCF header
samtools view -H "${CRAM}"  \
| grep '^@RG'  \
| awk '
{
  line=$0
  gsub(/\\/,"\\\\",line)   # escape backslashes first
  gsub(/\t/,"\\t",line)    # escape literal tabs
  print "##RG=" line
}' > readgroups.vcf.hdr


# Inject RG header lines into gvcf
# NOTE: Haplotypecaller ALWAYS outputs intervals in the GVCF, even if there are no reads - so drop these with bcftools
bcftools annotate \
  --header-lines readgroups.vcf.hdr \
  tmp.g.vcf.gz \
  | bcftools view \
    -e 'ALT="<NON_REF>" && (MAX(FORMAT/DP)=0 || MAX(FORMAT/MIN_DP)=0 || MAX(FORMAT/GQ)=0)' \
    -Oz -o "${IHASH}.${SAMPLE}.g.vcf.gz" 

bcftools index -t "${IHASH}.${SAMPLE}.g.vcf.gz" 

rm -f tmp.g.vcf.gz tmp.g.vcf.gz.tbi