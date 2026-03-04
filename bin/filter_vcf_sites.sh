#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed
# $6 = interval_hash
# $7 = pop

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# TODO: Break out vcf masks into individual components (i.e. Genmap, longdust, etc)

# Map variant type to bcftools selectors
# TODO: select variants on output of joint genotype
case "${4}" in
  snp)       TYPE_ARGS="-v snps   -m2 -M2 -e 'ALT=\"*\"'";;   # biallelic SNPs, drop star alleles
  indel)     TYPE_ARGS="-v indels -m2 -M2";;                  # biallelic INDELs
  invariant) TYPE_ARGS="-v ref";;                              # reference-only sites (if present)
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

# Subset to target variant class and run drop all genotypes, then run site-level soft filtering 
# (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)

# Subset to just the samples and update tags, then 
# Note MAC is calculated from MAF (7 decimal precision), this could cause rounding for very large cohorts (i.e. 100k+)
bcftools view --threads ${1} -S samples.list ${TYPE_ARGS} -Ou "${3}" \
  | bcftools +fill-tags -Ou - -- -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING,MAC=int(MAF*AN)' \
  | bcftools view -G -Ou \
  | bcftools filter -Ou -s QUAL_FAIL       -m+ -e "QUAL < ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s DP_FAIL         -m+ -e "INFO/DP < ${DPmin:-0} || INFO/DP < ${DPlower:-0} || INFO/DP > ${DPupper:-999999999}" \
  | bcftools filter -Ou -s DIST_INDEL_FAIL -m+ -e "INFO/DIST_INDEL < ${DIST_INDEL:--999999999}" \
  | bcftools filter -Ou -s EH_FAIL         -m+ -e "INFO/ExcHet < ${EH:--1}" \
  | bcftools filter -Ou -s HWE_FAIL        -m+ -e "INFO/HWE < ${HWE:--1}" \
  | bcftools filter -Ou -s MAF_FAIL        -m+ -e "INFO/MAF < ${MAF:-0}" \
  | bcftools filter -Ou -s MAC_FAIL        -m+ -e "INFO/MAC < ${MAC:-0}" \
  | bcftools filter -Ou -s NS_FAIL         -m+ -e "INFO/NS < ${NS:-0}" \
  | bcftools filter -Ou -s CR_FAIL         -m+ -e "INFO/CR < ${CR:-0}" \
  | bcftools filter -Ou -s MASK_FAIL       -m+ -M vcf_masks.bed \
  | bcftools view --threads ${1} -Oz -o ${4}.${7}.${6}.filt.vcf.gz
bcftools index --threads ${1} -t ${4}.${7}.${6}.filt.vcf.gz
