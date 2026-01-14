#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus
# $2 = mem (GB)
# $3 = outname

# Collect and sort files (initial list)
ls *.vcf.gz | sort -> vcf.list

# Detect vcf type (.g.vcf.gz or .vcf.gz) from the first file
first=$(head -n1 vcf.list || true)
if [[ -z "$first" ]]; then
    echo "No VCFs found (expected *.vcf.gz)" >&2
    exit 1
fi

# Exit if extension not recognised
if [[ "$first" == *.g.vcf.gz ]]; then
    ext=".g.vcf.gz"
elif [[ "$first" == *.vcf.gz ]]; then
    ext=".vcf.gz"
else
    echo "File extension not recognised: $first" >&2
    exit 1
fi

# Create list of vcfs in genomic order
# contig order from header
bcftools view -h "$first" \
| awk -F'[=,>]' '$0 ~ /^##contig=<ID=/ {print $3}' \
| awk '{print $1 "\t" NR}' \
> contig_rank.tsv

meta=$(mktemp)
skipped=$(mktemp)
: > "$meta"
: > "$skipped"

# build ordering keys from first record in each file; skip empties
while IFS= read -r f; do
  rec=$(bcftools view -H "$f" 2>/dev/null | head -n1 || true)

  # skip empty VCFs (no variant records)
  if [[ -z "$rec" ]]; then
    echo "$f" >> "$skipped"
    continue
  fi

  chrom=$(printf "%s\n" "$rec" | cut -f1)
  pos=$(printf "%s\n" "$rec" | cut -f2)

  rank=$(awk -v c="$chrom" '$1==c{print $2; found=1; exit} END{if(!found) print 999999}' contig_rank.tsv)
  printf "%s\t%s\t%s\n" "$rank" "$pos" "$f" >> "$meta"
done < vcf.list

# overwrite vcf.list with genomic-order list (empties removed)
LC_ALL=C sort -k1,1n -k2,2n "$meta" | cut -f3 > vcf.list

# Warn if some vcfs did not contain any records
if [[ -s "$skipped" ]]; then
  echo "Skipped empty VCFs: $(wc -l < "$skipped")" >&2
fi

# guard: everything empty, 
if [[ ! -s vcf.list ]]; then
  echo "All VCFs were empty; nothing to concat" >&2
  mv "${first}" "${3}${ext}" 
  bcftools index -t --threads ${1} "${3}${ext}" 
  exit 0
fi

# Run bcftools concat in --naive mode which is much faster when the chunks are already in global genomic order
# Use z9 for maximum compression, as these files will be stored long term
bcftools concat \
    --naive \
    --threads ${1} \
    -f vcf.list \
    -O z9 \
    -o "${3}${ext}"

# Index output - catch instances where naive concat has ruined contig order, in which case sort to fix
# NOTE This is triggered when re-merging indel,snp, and invariant vcfs
if bcftools index -t --threads ${1} "${3}${ext}" >/dev/null 2>&1; then
    echo "OK: VCF is sorted and indexable"
else
    echo "NOT OK: VCF is not sorted (or has other issues)"
    bcftools sort "${3}${ext}" --max-mem "${2}G" -Oz -o "${3}_sorted${ext}"
    mv -f "${3}_sorted${ext}" "${3}${ext}"
    bcftools index -t --threads ${1} "${3}${ext}"
fi
