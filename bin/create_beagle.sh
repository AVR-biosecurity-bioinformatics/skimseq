#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = gvcf
# $3 = ref_genome

# get interval hash
HASH=$( echo ${2%%.g.vcf.gz} )

# Remove <NON_REF> alleles at sites with variants but retain them at HOM REF
# then remove multiallelic sites, indel sites, sites with no depth, and spanning indels '*'
# NOTE should i replace sites with star alleles to <NON_REF> rather than filtering them out altogether?
bcftools view -A $2 \
    | bcftools view \
        -m2 \
        -M2 \
        --exclude-types indels \
        -e 'sum(FMT/DP) == 0 || ALT=="*"' \
        -o filtered_output.vcf

# Count the types of variants that made it through
#bcftools query -f '%ALT\n' filtered_output.vcf | sort | uniq -c
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' filtered_output.vcf 

# PP fields isnt filled for missing data, so replace with a uniform prior
bcftools query -f '%CHROM\t%POS[\t%PP]' filtered_output.vcf \
    | awk -F'\t' 'BEGIN {OFS="\t"} {for(i=3; i<=NF; i++) if ($i == ".") $i="0,0,0"; print $0}' \
    | bgzip \
    > updated_pp.txt.gz

# Index the file with tabix
tabix -s1 -b2 -e2 updated_pp.txt.gz

# Replace the PP field with the updated ones
bcftools annotate \
    -a updated_pp.txt.gz \
    -c CHROM,POS,FORMAT/PP \
    filtered_output.vcf \
    -o filtered_output_updated.vcf

# Drop existing PL tag and change the name of the PP tag to PL 
# NOTE - Skip this step if a genotype likelihood rather than posterior file is required
echo FORMAT/PP PL > rename_file
bcftools annotate \
    -x FORMAT/PL \
    --rename-annots rename_file \
    filtered_output_updated.vcf \
    > filtered_output_posteriors.vcf

# Use BCFtools to convert phred scaled likelihoods into probabilities (similar to beagle file)
bcftools +tag2tag \
    filtered_output_posteriors.vcf \
    -- \
    -r \
    --pl-to-gp \
    > filtered_rescaled.vcf

# Expand gvcf reference blocks
bcftools convert \
    --gvcf2vcf filtered_rescaled.vcf \
    --fasta-ref $3 \
    > expanded.vcf

# Create beagle file of posterior probabilities
paste -d '\t' \
    <(echo -e "marker\tallele1\tallele2") \
    <(bcftools query -l expanded.vcf \
    | awk '{for(i=1; i<=3; i++) print $0}' \
    | tr '\n' '\t') \
    > ${HASH}.beagle

bcftools query \
    -f '%CHROM:%POS\t%REF\t%ALT[\t%GP]\n' \
    expanded.vcf \
    | tr ',' '\t' \
    >> ${HASH}.beagle