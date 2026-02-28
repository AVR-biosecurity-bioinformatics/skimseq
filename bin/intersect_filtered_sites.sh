#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = global_vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = interval_hash

# Subset to target variant class and run drop all genotypes, then run site-level soft filtering 
# (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)

# threshold: exclude if failed in > N pops
N=2   # TODO: make this a nextflow param passed in

# 1) build BED of positions failing in >N pops
tmp_fail_bed=$(mktemp)

# count failures by CHROM,POS across pop vcfs
# (requires pop vcfs to include FILTER field with FAIL tags; not PASS-only)
while read -r v; do
  [[ -z "$v" ]] && continue
  bcftools query -i 'FILTER!="PASS"' -f '%CHROM\t%POS\n' "$v"
done < pop_vcfs.list \
| LC_ALL=C sort -k1,1 -k2,2n \
| uniq -c \
| awk -v n="$N" 'BEGIN{OFS="\t"} $1>n {chr=$2; pos=$3; print chr, pos-1, pos}' \
> "$tmp_fail_bed"

# 2) filter the global sitelist: keep PASS sites (if global is tagged), then exclude fail bed
out="${4}.${6}.final.sites.vcf.gz"

bcftools view --threads "${1}" -f PASS -Ou "${3}" \
  | bcftools view --threads "${1}" -T ^"$tmp_fail_bed" -Oz -o "$out"

bcftools index --threads "${1}" -t "$out"

# counts
nvars=$(bcftools index -n "$out" | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${6}.counts"

# summary: how many sites were excluded for failing in >N pops
excluded=$(awk 'END{print NR+0}' "$tmp_fail_bed")
printf "EXCLUDED_FAIL_GT%d_POPS\t%d\n" "$N" "$excluded" > "${4}.${6}_filter_summary.tsv"

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n "${4}.${6}.sites.vcf.gz" | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${6}.counts"

# Create a small summary of the number of sites passing and failing each filter
bcftools query -f '%FILTER\n' $tmp_fail_bed \
  | sort \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"} {print $2, $1}' \
  > "${4}.${6}_filter_summary.tsv"

# Remove temporary vcf files
rm -f "$tmp_fail_bed"
