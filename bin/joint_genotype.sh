#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = genomicsDB
# $4 = ref_genome
# $5 = interval hash
# $6 = interval_bed
# $7 = exclude_bed
# $8 = exclude_padding
# $9 = output_invariant

# joint genotype variant sites only
gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
    -R ${4} \
    -V gendb://${3} \
    -L ${6} \
    -O ${5}.vcf.gz \
    --exclude-intervals ${7} \
    --interval-exclusion-padding ${8} \
    --interval-merging-rule ALL \
    --merge-input-intervals true \
    --only-output-calls-starting-in-intervals \
    --include-non-variant-sites ${9} \
    --tmp-dir /tmp

#    --use-posteriors-to-calculate-qual true \
#    --genotype-assignment-method USE_POSTERIORS_ANNOTATION \