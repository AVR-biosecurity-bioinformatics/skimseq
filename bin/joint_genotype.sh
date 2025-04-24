#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = gvcf
# $3 = ref_genome
# $4 = interval_list

# get interval hash
HASH=$( echo ${4%%.interval_list} )

# joint genotype variant sites only
gatk --java-options "-Xmx8G" GenotypeGVCFs \
    -R $3 \
    -V $2 \
    -O ${HASH}.vcf.gz \
    --use-posteriors-to-calculate-qual true \
    --genotype-assignment-method USE_POSTERIORS_ANNOTATION \
    --tmp-dir /tmp