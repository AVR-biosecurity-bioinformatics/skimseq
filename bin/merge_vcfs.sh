#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus
# $2 = mem (GB)
# $3 = outname

# Collect and sort files in the new numbered format (this assumes files are like _00001.vcf.gz or _00001.g.vcf.gz)
ls _*.vcf.gz | sort -V > vcf.list

# 2. detect type from the first file
first=$(head -n1 vcf.list || true)
if [[ -z "$first" ]]; then
    echo "No VCFs found (expected _*.vcf.gz)" >&2
    exit 1
fi

if [[ "$first" == *.g.vcf.gz ]]; then
    ext=".g.vcf.gz"
elif [[ "$first" == *.vcf.gz ]]; then
    ext=".vcf.gz"
else
    echo "File extension not recognised: $first" >&2
    exit 1
fi

# Run bcftools concat in naive mode, --naive because the chunks are already in global genomic order
# Use z9 for maximum compression, as these files will be stored long term
bcftools concat \
    --naive \
    --threads ${1} \
    -f vcf.list \
    -O z9 \
    -o "${3}${ext}"

# index the output (tabix)
bcftools index -t --threads "$CPUS" "${3}${ext}"