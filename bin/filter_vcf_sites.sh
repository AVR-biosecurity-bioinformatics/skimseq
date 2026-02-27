#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed
# $6 = interval_hash
# $7 = DP summary

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

# Calculate percentile DP filters from DP histogram
read DPlower DPupper < <(
  awk -v pl="$PCT_LOW" -v ph="$PCT_HIGH" '
    { dp[NR]=$1; cnt[NR]=$2+0; N+=cnt[NR] }
    END{
      if(N==0) exit 1
      low  = pl/100*N; li=int(low);  if(li<low)  li++; if(li<1) li=1; if(li>N) li=N
      high = ph/100*N; ui=int(high); if(ui<high) ui++; if(ui<1) ui=1; if(ui>N) ui=N
      cum=0
      for(i=1;i<=NR;i++){
        cum += cnt[i]
        if(!lo && cum>=li) lo=dp[i]
        if(!hi && cum>=ui){ hi=dp[i]; break }
      }
      printf "%d %d\n", lo, hi
    }
  ' "$7"
)

# Subset to target variant class and run drop all genotypes, then run site-level soft filtering 
# (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)

# TODO: subset to samples here to allow for pop-specific calculation
bcftools view --threads ${1} -G ${TYPE_ARGS} -Ou "${3}" \
  | bcftools filter -Ou -s QUAL_FAIL       -m+ -e "QUAL <= ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s DP_FAIL         -m+ -e "INFO/DP <= ${DPmin:-0} || INFO/DP <= ${DPlower:-0} || INFO/DP >= ${DPupper:-999999999}" \
  | bcftools filter -Ou -s DIST_INDEL_FAIL -m+ -e "INFO/DIST_INDEL <= ${DIST_INDEL:--999999999}" \
  | bcftools filter -Ou -s EH_FAIL         -m+ -e "INFO/ExcHet <= ${EH:-1e9}" \
  | bcftools filter -Ou -s HWE_FAIL        -m+ -e "INFO/HWE <= ${HWE:-1e9}" \
  | bcftools filter -Ou -s MAF_FAIL        -m+ -e "INFO/MAF <= ${MAF:-0}" \
  | bcftools filter -Ou -s MAC_FAIL        -m+ -e "INFO/MAC <= ${MAC:-0}" \
  | bcftools filter -Ou -s NS_FAIL         -m+ -e "INFO/NS <= ${NS:-0}" \
  | bcftools filter -Ou -s CR_FAIL         -m+ -e "INFO/CR <= ${CR:-0}" \
  | bcftools filter -Ou -s MASK_FAIL       -m+ -M vcf_masks.bed \
  | bcftools view --threads ${1} -Ob -o tmp.tagged.bcf

# Keep only variants that PASS & index output
# TODO: Drop FT and other extra fields from vcf
bcftools view --threads ${1} -f PASS -Oz -o ${4}.${6}.sites.vcf.gz tmp.tagged.bcf
bcftools index --threads ${1} -t ${4}.${6}.sites.vcf.gz

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n "${4}.${6}.sites.vcf.gz" | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${6}.counts"

# Create a small summary of the number of sites passing and failing each filter
bcftools query -f '%FILTER\n' tmp.tagged.bcf \
  | sort \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"} {print $2, $1}' \
  > "${4}.${6}_filter_summary.tsv"

# Remove temporary vcf files
rm -f tmp* *.bcf