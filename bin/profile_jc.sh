#!/bin/bash
set -euo pipefail


# Args
cpus="${1:-1}"
mem="${2:-0}"                # placeholder
ihash="${3}"                 # interval_hash
logfile="${4}"               # GenotypeGVCFs stderr
vcf="${5}"                   # joint-called VCF/BCF (bgzipped & indexed)

tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

###############################################################################
# 1) Parse ProgressMeter → progress.tsv
###############################################################################
# Header: interval_hash chrom pos elapsed_min variants variants_per_min
printf "interval_hash\tchrom\tpos\telapsed_min\tvariants\tvariants_per_min\n" > "$tmpdir/progress.tsv"

grep -F 'ProgressMeter -' "$logfile" \
| awk -v h="$ihash" '
  # ... ProgressMeter -   CM028320.1:49552   0.0   24000   5347.7
  match($0,/ProgressMeter -[[:space:]]*([[:alnum:]_.]+):([0-9]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9.]+)/, m){
    printf("%s\t%s\t%s\t%.3f\t%d\t%.1f\n", h, m[1], m[2], m[3], m[4], m[5])
  }' >> "$tmpdir/progress.tsv"

###############################################################################
# 2) Build tick windows per (interval_hash, chrom)
#    First window: [0, first_pos]
###############################################################################
awk -F'\t' -v OFS='\t' '
NR==1{next}
{
  h=$1; c=$2; p=$3+0; em=$4+0; v=$5+0;
  key=h "|" c
  if(!(key in seen)){
    ws=0; we=(p>0?p:1); dem=em; dv=v; idx[key]=1
    print c, ws, we, h, dem, dv, (h "|" c "|" idx[key])
  } else {
    ws=prev_p[key]-1; if(ws<0) ws=0
    we=(p>ws?p:ws+1)
    dem=em-prev_em[key]; if(dem<0) dem=0
    dv =v -prev_v[key];  if(dv <0) dv =0
    idx[key]++
    print c, ws, we, h, dem, dv, (h "|" c "|" idx[key])
  }
  seen[key]=1; prev_p[key]=p; prev_em[key]=em; prev_v[key]=v
}' "$tmpdir/progress.tsv" > "$tmpdir/tick_windows.bed"
# cols: 1 chrom 2 start 3 end 4 interval_hash 5 d_elapsed_min 6 d_variants 7 window_id

###############################################################################
# 3) Make per-site annotations from VCF (one pass each)
#    - ALL sites (point BED)
#    - SNP sites
#    - INDEL sites
#    - site-level: n_alleles, NS (samples with data)
###############################################################################
# All sites
bcftools view -H "$vcf" \
| awk -v OFS='\t' '{print $1, $2-1, $2}' > "$tmpdir/sites_all.bed"

# SNPs
bcftools view -H -v snps "$vcf" \
| awk -v OFS='\t' '{print $1, $2-1, $2}' > "$tmpdir/sites_snp.bed"

# INDELs
bcftools view -H -v indels "$vcf" \
| awk -v OFS='\t' '{print $1, $2-1, $2}' > "$tmpdir/sites_indel.bed"

# Ensure INFO/NS exists; then extract per-site NS and number of alleles
# sum_alleles per site = 1 (REF) + #ALTs (count commas in ALT; ALT=="." → 1)
bcftools +fill-tags -Ou -t NS "$vcf" \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/NS\n' \
| awk -v OFS='\t' '
  {
    chrom=$1; pos=$2; ref=$3; alt=$4; ns=$5+0;
    if (alt=="." || alt=="") nalt=0;
    else { nalt=gsub(/,/, "&", alt) + 1 }   # commas+1 ALTs
    nalleles = 1 + nalt;                    # REF + ALTs
    print chrom, pos-1, pos, nalleles, ns
  }' > "$tmpdir/site_alleles_ns.bed"
# cols: chrom start end n_alleles NS

###############################################################################
# 4) Attribute sites to windows & summarise per window
###############################################################################
# Count variants/snp/indel per window (should match Δvariants)
bedtools intersect -a "$tmpdir/tick_windows.bed" -b "$tmpdir/sites_all.bed"   -c > "$tmpdir/_w_all.tsv"
bedtools intersect -a "$tmpdir/_w_all.tsv"        -b "$tmpdir/sites_snp.bed" -c > "$tmpdir/_w_snp.tsv"
bedtools intersect -a "$tmpdir/_w_snp.tsv"        -b "$tmpdir/sites_indel.bed" -c > "$tmpdir/_w_counts.tsv"
# _w_counts.tsv cols now:
# 1 chrom 2 start 3 end 4 ihash 5 d_elapsed 6 d_variants(pm) 7 window_id 8 n_variants(vcf) 9 n_snps 10 n_indels

# Sum n_alleles and NS per window
bedtools map -a "$tmpdir/tick_windows.bed" -b "$tmpdir/site_alleles_ns.bed" -c 4,5 -o sum,sum \
> "$tmpdir/_w_alleles_ns.tsv"
# cols: chrom start end ihash d_elapsed d_variants window_id sum_n_alleles sum_NS

# Join counts + alleles/ns on (window_id)
awk -F'\t' 'NR==FNR{a[$7]=$0; next}
{
  key=$7
  if(key in a){
    print a[key] "\t" $8 "\t" $9
  }
}' "$tmpdir/_w_counts.tsv" "$tmpdir/_w_alleles_ns.tsv" \
| awk -F'\t' -v OFS='\t' '
BEGIN{
  print "interval_hash","window_id","chrom","start","end","window_elapsed_min",
        "variants_pm","n_variants","n_snps","n_indels",
        "sum_alleles","genotypes_with_data"
}
{
  # From joined row:
  # counts: 1 chrom 2 start 3 end 4 ihash 5 d_elapsed 6 d_variants_PM 7 win_id 8 n_variants 9 n_snps 10 n_indels
  # alleles/ns: appended as 11 sum_n_alleles 12 sum_NS
  print $4,$7,$1,$2,$3,$5, $6,$8,$9,$10, $11,$12
}' > {3}.profile.tsv