#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = ref_genome
# $4 = posterior

# Note - beagle output only works for polymorphic sites but indels and nonvariants are supported

if [[ ${4} == "true" ]]; then
    # If posterior is true, replace PL tag with the PP tag for compatibility with bcftools
    echo FORMAT/PP PL > rename_file
    bcftools annotate \
        -x FORMAT/PL \
        --rename-annots rename_file \
        ${2} \
        -o tmp.vcf
    outname='pp'
else
    bcftools view ${2} \
        -o tmp.vcf
    outname='gl'
fi

# Use BCFtools to convert phred scaled likelihoods into probabilities (similar to beagle file)
bcftools +tag2tag \
    tmp.vcf  \
    -- \
    -r \
    --pl-to-gp \
    > tmp_rescaled.vcf

# Create beagle file header
paste -d '\t' \
    <(echo -e "marker\tallele1\tallele2") \
    <(bcftools query -l tmp_rescaled.vcf \
    | awk '{for(i=1; i<=3; i++) print $0}' \
    | tr '\n' '\t') \
    > ${outname}.beagle

# Fill in beagle file
bcftools query \
    -f '%CHROM:%POS\t%REF\t%ALT[\t%GP]\n' \
    tmp_rescaled.vcf \
    | tr ',' '\t' \
    >> ${outname}.beagle

rm -f tmp*.vcf*