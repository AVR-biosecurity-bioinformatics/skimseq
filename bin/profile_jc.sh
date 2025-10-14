#!/bin/bash
set -euo pipefail

## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = one or more *.stderr.log files

# Collect files exactly as given after $2
LOGS=( "${@:3}" )

PROG_OUT="progress_summary.tsv"
ALLELE_OUT="too_many_alleles.tsv"
SKIP_OUT="skipped_sites.tsv"

# headers
printf "source\tchrom\tpos\telapsed_min\tvariants\tvariants_per_min\n" > "$PROG_OUT"
printf "source\tchrom\tpos\n" > "$ALLELE_OUT"
printf "source\tchrom\tpos\n" > "$SKIP_OUT"

for log in "${LOGS[@]}"; do
  # Source label: basename, drop leading underscore (optional), drop suffix
  base="$(basename -- "$log")"
  src="${base#_}"              # remove leading '_' if present; delete this line to keep it
  src="${src%.stderr.log}"

  # Progress lines → merged TSV
  grep -F 'ProgressMeter -' "$log" 2>/dev/null | \
  awk -v src="$src" '
    /ProgressMeter -/ {
      chrompos=$3; elapsed=$4; nvar=$5; rate=$6
      split(chrompos,a,":"); chrom=a[1]; pos=a[2]
      if (elapsed ~ /^[0-9.]+$/ && nvar ~ /^[0-9]+$/ && rate ~ /^[0-9.]+$/)
        printf("%s\t%s\t%s\t%.3f\t%d\t%.1f\n", src, chrom, pos, elapsed, nvar, rate)
    }' >> "$PROG_OUT" || true

  # Too-many-alleles positions → merged TSV
  grep -F 'has too many alleles in the combined VCF record' "$log" 2>/dev/null | \
  awk -v src="$src" '{
    chrom=""; pos=""
    for(i=1;i<=NF;i++){
      if($i=="Chromosome"){chrom=$(i+1)}
      if($i=="position"){pos=$(i+1)}
    }
    sub(/:$/,"",chrom); gsub(/[^0-9]/,"",pos)
    if(chrom!="" && pos!="") printf("%s\t%s\t%s\n", src, chrom, pos)
  }' >> "$ALLELE_OUT" || true

  # Skipped sites (insufficient data) → merged TSV
  grep -F 'MinimalGenotypingEngine - Some genotypes contained insufficient data' "$log" 2>/dev/null | \
  awk -v src="$src" -F'location ' '
    NF>1 {
      split($2,a,":"); chrom=a[1]; pos=a[2]; gsub(/[[:space:]]/,"",pos)
      if(chrom!="" && pos!="") printf("%s\t%s\t%s\n", src, chrom, pos)
    }' >> "$SKIP_OUT" || true
done

printf "Wrote: %s, %s, %s\n" "$PROG_OUT" "$ALLELE_OUT" "$SKIP_OUT"