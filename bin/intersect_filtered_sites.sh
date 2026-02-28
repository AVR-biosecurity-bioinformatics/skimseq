#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = global_vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = interval_hash
# $6 = n_pops_failing
# $7 = perc_pops_failing

# Subset to target variant class and run drop all genotypes, then run site-level soft filtering 
# (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)

# number of pops (one vcf per pop)
n_pops=$(awk 'NF{n++} END{print n+0}' pop_vcfs.list)

# percentage-based threshold: Tp = ceil(PERC * n_pops / 100)
# (use ceil so 1% of 3 pops => 1, not 0; change to floor if you prefer)
perc_N=$(awk -v p="${7}" -v n="$n_pops" '
  BEGIN{
    if(p<=0 || n<=0){print 0; exit}
    x = p*n/100.0
    t = int(x)
    if(t < x) t++
    print t
  }'
)

# effective threshold is the larger of N and percent-based
N=${6}
if (( perc_N > N )); then N=$perc_N; fi

# 1) build BED of positions failing in >N pops
tmp_fail_bed=$(mktemp)

# count failures by CHROM,POS across pop vcfs
while read -r v; do
  [[ -z "$v" ]] && continue
  bcftools query -i 'FILTER!="PASS"' -f '%CHROM\t%POS\n' "$v"
done < pop_vcfs.list \
| LC_ALL=C sort -k1,1 -k2,2n \
| uniq -c \
| awk -v n="$N" 'BEGIN{OFS="\t"} $1>n {chr=$2; pos=$3; print chr, pos-1, pos}' \
> "$tmp_fail_bed"

# 2) filter the global sitelist: keep PASS sites (if global is tagged), then exclude fail bed
out="${4}.${5}.sites.vcf.gz"

bcftools view --threads "${1}" -f PASS -Ou "${3}" \
  | bcftools view --threads "${1}" -T ^"$tmp_fail_bed" -Oz -o "$out"

bcftools index --threads "${1}" -t "$out"

# summary: how many sites were excluded for failing in >N pops
excluded=$(awk 'END{print NR+0}' "$tmp_fail_bed")
printf "EXCLUDED_FAIL_GT%d_POPS\t%d\n" "$N" "$excluded" > "${4}.${5}_filter_summary.tsv"

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n $out | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${5}.counts"

# Remove temporary files
#rm -f "$tmp_fail_bed"
