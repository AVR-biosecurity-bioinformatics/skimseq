#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# Map VTYPE -> bcftools selectors
case "${4}" in
  snp)       TYPE_ARGS="-v snps   -m2 -M2 -e 'ALT=\"*\"'";;   # biallelic SNPs, drop star alleles
  indel)     TYPE_ARGS="-v indels -m2 -M2";;                  # biallelic INDELs
  invariant) TYPE_ARGS="-v ref";;                              # reference-only sites (if present)
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

## TODO: Move all annotations to here

# Subset to target variant class
bcftools view --threads ${1} ${TYPE_ARGS} -Ou "${3}" \
| \
# Filter genotypes -> set failing GTs to missing
bcftools +setGT --threads ${1} -Ou -- -t q -n . \
    -i "FMT/GQ < ${GQ:-0} || FMT/DP < ${gtDPmin:-0} || FMT/DP > ${gtDPmax:-999999}" \
| \
# Recompute site tags that depend on GTs
bcftools +fill-tags --threads ${1} -Ou -- -t AC,AN,AF,MAF \
| \
# Site-level filtering (uses env vars exported by Nextflow)
bcftools filter --threads ${1} -Ou -e "
    (INFO/QD < ${QD:-0}) ||
    (QUAL < ${QUAL_THR:-0}) ||
    (INFO/SOR > ${SOR:-1e9}) ||
    (INFO/FS  > ${FS:-1e9})  ||
    (INFO/MQ  < ${MQ:-0})    ||
    (INFO/MQRankSum      < ${MQRS:--1e9}) ||
    (INFO/ReadPosRankSum < ${RPRS:--1e9}) ||
    (INFO/MAF < ${MAF:-0}) ||
    (INFO/MAC < ${MAC:-0}) ||
    (INFO/ExcessHet > ${EH:-1e9}) ||
    (INFO/DP < ${DPmin:-0}) ||
    (INFO/DP > ${DPmax:-999999999}) ||
    (F_MISSING > ${F_MISSING:-1})
" \
-M vcf_masks.bed \
--soft-filter + \
| \
# Write output + index
bcftools view --threads ${1} -Oz -o ${4}_filtered.vcf.gz
bcftools index --threads ${1} -t ${4}_filtered.vcf.gz

# TODO: Add exit code passing for piper