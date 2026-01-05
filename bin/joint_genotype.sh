#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = genomicsDB
# $4 = ref_genome
# $5 = interval hash
# $6 = interval_bed
# $7 = exclude_bed

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
GENOMICSDB="${3}"
REF="${4}"
IHASH="${5}"
INTERVAL_BED="${6}"
EXCLUDE_BED="${7}"

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${MEM_GB} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
# Send stderr to log file for profiling
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" GenotypeGVCFs \
    -R "${REF}" \
    -V gendb://${GENOMICSDB} \
    -L "${INTERVAL_BED}" \
    -O ${IHASH}.vcf.gz \
    --exclude-intervals "${EXCLUDE_BED}" \
    --interval-exclusion-padding "${EXCLUDE_PAD}" \
    --include-non-variant-sites "${OUTPUT_INVARIANT}" \
    --interval-merging-rule ALL \
    --merge-input-intervals \
    --variant-output-filtering STARTS_IN \
    --max-alternate-alleles "${MAX_ALTERNATE}" \
    --genomicsdb-max-alternate-alleles "${GENOMICSDB_MAX_ALTERNATE}" \
    -ploidy "${PLOIDY}" \
    --heterozygosity "${HET}" \
    --heterozygosity-stdev "${HET_SD}" \
    --indel-heterozygosity "${INDEL_HET}" \
    --tmp-dir /tmp \
    --genomicsdb-shared-posixfs-optimizations true \
    2> >(tee -a ${IHASH}.stderr.log >&2)