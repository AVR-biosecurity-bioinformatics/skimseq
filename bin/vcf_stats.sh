#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = Reference genome
# $4 = Sample

# Subset to sample and exclude uncalled
bcftools view -s ${4} --exclude-uncalled -Ob -o ${4}.bcf ${2}
bcftools index ${4}.bcf

# Calculate stats
bcftools stats -F ${3} ${4}.bcf \
| awk -v s="${4}" 'BEGIN{FS=OFS="\t"}
    /^#/ {print; next}
    $1=="ID" { $3=s ".vcf.gz" }
    { print }
' > "${4}.vcfstats.txt"

rm ${4}.bcf*