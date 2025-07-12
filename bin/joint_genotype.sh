#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = genomicsDB
# $4 = ref_genome
# $5 = interval hash
# $6 = interval_list

# joint genotype variant sites only
gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
    -R ${4} \
    -V gendb://${3} \
    -L ${6} \
    -O ${5}.vcf.gz \
    --tmp-dir /tmp

#    --use-posteriors-to-calculate-qual true \
#    --genotype-assignment-method USE_POSTERIORS_ANNOTATION \