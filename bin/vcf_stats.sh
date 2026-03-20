#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = Reference genome
# $4 = Sample

bcftools stats \
    --threads ${1} \
    -F ${3} \
    -s - \
    ${2} > "vcfstats.txt"

#bcftools view \
#  --threads ${1} \
#  -s ${4} \
#  --exclude-uncalled \
#  -Ou ${2} \
#| bcftools stats \
#    --threads ${1} \
#    -F ${3} \
#    -s ${4} \
#    - \
#| awk -v s="${4}" 'BEGIN{FS=OFS="\t"}
#    /^#/ {print; next}
#    $1=="ID" { $3=s ".vcf.gz" }
#    { print }
#' > "${4}.vcfstats.txt"