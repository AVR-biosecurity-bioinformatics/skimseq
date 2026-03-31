#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus
# $2 = mem (GB)
# $3 = vcf
# $4 = outname

 bcftools view -Ou ${3} \
| tee \
    >(bcftools view -Oz -v snps   -o ${4}.snp.vcf.gz) \
    >(bcftools view -Oz -v indels -o ${4}.indel.vcf.gz) \
| bcftools view -Oz -v ref -o ${4}.invariant.vcf.gz

# Index outputs
bcftools index -t ${4}.snp.vcf.gz
bcftools index -t ${4}.indel.vcf.gz
bcftools index -t ${4}.invariant.vcf.gz
