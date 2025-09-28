#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# Map variant type to bcftools selectors
# TODO: select variants on output of joint genotype
case "${4}" in
  snp)       TYPE_ARGS="-v snps   -m2 -M2 -e 'ALT=\"*\"'";;   # biallelic SNPs, drop star alleles
  indel)     TYPE_ARGS="-v indels -m2 -M2";;                  # biallelic INDELs
  invariant) TYPE_ARGS="-v ref";;                              # reference-only sites (if present)
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

# Subset to target variant class
bcftools view --threads ${1} ${TYPE_ARGS} -Ou "${3}" \
| \
# Filter genotypes -> set failing GTs to missing
bcftools +setGT --threads ${1} -Ou -- -t q -n . \
    -i "FMT/GQ < ${GQ:-0} || FMT/DP < ${gtDPmin:-0} || FMT/DP > ${gtDPmax:-999999}" \
| \
# TODO: filter samples using -S samples_to_keep.txt
bcftools view --threads ${1} -U -o tmp.vcf.gz 

# Recompute site tags that depend on GTs and samples
bcftools +fill-tags --threads ${1} tmp.vcf.gz -o tmp2.vcf.gz -- -t AC,AN,MAF,F_MISSING,NS

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp2.vcf.gz \
| awk 'BEGIN{OFS="\t"}
       {
         split($3,ac,",")          # AC is comma‑separated if multi‑allelic
         mac=$4                    # start with AN
         refCount = $4             # will be AN - sum(AC)
         for(i in ac){refCount-=ac[i]; mac=(ac[i]<mac?ac[i]:mac)}
         mac=(refCount<mac?refCount:mac)
         print $1,$2,mac
       }'                         \
| bgzip > MAC.tsv.gz
tabix -s1 -b2 -e2 MAC.tsv.gz

# Create VCF header line for MAC filter
cat > add_MAC.hdr <<'EOF'
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count (minimum of each ALT AC and reference allele count)">
EOF

# Annotate the vcf with INFO/MAC
bcftools annotate \
     -h add_MAC.hdr \
     -a MAC.tsv.gz \
     -c CHROM,POS,INFO/MAC \
     -O z -o tmp3.vcf.gz \
     tmp2.vcf.gz

# Site-level filtering (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)
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
--soft-filter FILTER \
tmp3.vcf.gz \
| \
# Keep only variants that PASS
bcftools view --threads ${1} -f PASS -Oz -o ${4}_filtered.vcf.gz

# Index output vcf
bcftools index --threads ${1} -t ${4}_filtered.vcf.gz

# Remove temporary vcf files
rm -f tmp* *.tsv.gz