#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref_genome
# $4 = interval hash
# $5 = targets_file

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
IHASH="${4}"
TARGETS_FILE="${5}"

# Set up intervals - handle fixed sitelists

if [[ "${TARGETS_FILE}" == *.vcf.gz ]]; then
  echo "[targets] Detected VCF panel: ${TARGETS_FILE}" >&2

  ## Build allele targets: CHROM POS REF,ALT  (tabix indexed)
  bcftools view -m2 -M2 "${TARGETS_FILE}" \
    | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' \
    | bgzip -c > panel.alleles.tsv.gz
  tabix -s1 -b2 -e2 panel.alleles.tsv.gz

  # Build a BED (0-based) for filtering CRAMs by region
  bcftools view -m2 -M2 "${TARGETS_FILE}" \
    | bcftools query -f'%CHROM\t%POS0\t%POS\n' \
    | bgzip -c > panel.bed.gz
  tabix -p bed panel.bed.gz

  INTERVAL_BED="panel.bed.gz"

  # Use allele-restricted calling at those positions and ensure all are present
  MPILEUP_TARGETS_FLAGS="-T panel.alleles.tsv.gz"
  CALL_FLAGS="-C alleles -T panel.alleles.tsv.gz --insert-missed"

elif [[ "${TARGETS_FILE}" == *.bed* ]]; then
  echo "[targets] Detected BED intervals: ${TARGETS_FILE}" >&2
  # For BED, just target the intervals; no allele restriction possible
  MPILEUP_TARGETS_FLAGS="-R ${TARGETS_FILE}"
  INTERVAL_BED=${TARGETS_FILE}
  CALL_FLAGS=""
else
  echo "${TARGETS_FILE} is in wrong format" >&2
  exit 1
fi 

# -----------------------------
# Pre-filter reads using samtools
# -----------------------------

# Set up filtering expressions for samtools

# Min aligned bases (excluding softclips)
ALN_CLAUSE="(qlen - sclen) >= ${MIN_ALIGNED_LENGTH}"

# Min and max fragment (insert size) length
TLEN_CLAUSE="( (tlen >= ${MIN_FRAGMENT_LENGTH} && tlen <= ${MAX_FRAGMENT_LENGTH}) || (tlen <= -${MIN_FRAGMENT_LENGTH} && tlen >= -${MAX_FRAGMENT_LENGTH}) )"

# Read flags filter (make DUP conditional)
if [[ "${RMDUP}" == "false" ]]; then
  FLAGS_CLAUSE='(!flag.unmap && !flag.secondary && !flag.supplementary  && flag.proper_pair )'
else
  FLAGS_CLAUSE='(!flag.unmap && !flag.secondary && !flag.supplementary  && flag.proper_pair && !flag.dup)'
fi

# Create joint filtering expression for samtools
EXPR="${FLAGS_CLAUSE} && ${ALN_CLAUSE} && ${TLEN_CLAUSE}"

export REF INTERVAL_BED EXPR

# Choose jobs and cpus for GNU parallel
THREADS_PER_JOB=2
JOBS=$(( CPUS / THREADS_PER_JOB )) # how many CRAMs in parallel
(( JOBS < 1 )) && JOBS=1

# Subset and filter CRAMs in parallel
parallel --jobs "${JOBS}" --line-buffer '
  cram={}
  out="${cram%.cram}.filt.cram"

  samtools view -@ '"${THREADS_PER_JOB}"' -T "${REF}" \
    --regions-file "${INTERVAL_BED}" \
    -e "${EXPR}" \
    -O cram -o "${out}" "${cram}"

  samtools index -@ '"${THREADS_PER_JOB}"' "${out}"
  echo "${out}"
' :::: cram.list | sort > cram.filtered.list

# -----------------------------
# Variant calling using pre-filtered crams
# -----------------------------

# set up bcftools filter flags
if [[ "${RMDUP}" == "false" ]]; then
  FILTER_FLAGS="--ns DUP -G UNMAP,SECONDARY,QCFAIL"
else
  FILTER_FLAGS="--ns UNMAP,SECONDARY,QCFAIL,DUP"
fi

if [[ "${OUTPUT_INVARIANT}" == "false" ]]; then
  VARIANTS_ONLY="--variants-only"
else
  VARIANTS_ONLY=""
fi

# Call variants with mpleup
bcftools mpileup \
    --threads ${CPUS} \
    --bam-list cram.filtered.list \
    --max-depth ${MAXDEPTH} \
    --fasta-ref ${REF} \
    --min-BQ ${MINBQ} \
    --min-MQ ${MINMQ} \
    ${MPILEUP_TARGETS_FLAGS} \
    ${FILTER_FLAGS} \
    --annotate FORMAT/DP,FORMAT/AD,INFO/AD \
    --indels-cns \
    --indel-size 110 \
    | bcftools call \
    -Ou \
    -a FORMAT/GP,FORMAT/GQ \
    --ploidy ${PLOIDY} \
    ${CALL_FLAGS} \
    ${VARIANTS_ONLY} \
    --multiallelic-caller \
    --prior ${MUTATION_RATE} \
  | bcftools +setGT \
    -Ou -- \
    -t q -n . -i 'FMT/DP=0' \
  | bcftools annotate \
    --threads "${CPUS}" \
    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    -Oz9 -o "${IHASH}.vcf.gz"

# index output
bcftools index -t ${IHASH}.vcf.gz

# Extra annotatiomns that can be added
# mpileup: FORMAT/SP
# call: FORMAT/GP,INFO/PV4

# Clean up temporary files
xargs -r -d '\n' rm -f < <(awk '{print $0; print $0 ".crai"}' cram.filtered.list)
rm -f cram.filtered.list