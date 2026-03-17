#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome
# $5 = Sample

bcftools view \
  --threads ${1} \
  -s ${5} \
  --exclude-uncalled \
  -Ou ${2} \
| bcftools stats \
    --threads ${1} \
    -F ${4} \
    -s ${5} \
    - \
| awk -v s="${5}" 'BEGIN{FS=OFS="\t"}
    /^#/ {print; next}
    $1=="ID" { $3=s ".vcf.gz" }
    { print }
' > "${5}.vcfstats.txt"