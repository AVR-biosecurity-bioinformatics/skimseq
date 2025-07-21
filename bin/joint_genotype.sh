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

# First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
if [[${9} == "false" ]]; then
    # joint genotype variant sites only
    gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O joint_called.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals true \
        --only-output-calls-starting-in-intervals \
        --include-non-variant-sites false \
        --tmp-dir /tmp
else 
    # Joint genotype both variant and invariant
    # This requires some custom code to re-add missing genotpye fields for compatibility with later steps

    # First use gatk selectvariants to get the sites to genotype
    # This resolves the memory leak when calling invariant sites from genomicsDB reported in: https://github.com/broadinstitute/gatk/issues/8989
    gatk --java-options "-Xmx${2}G" SelectVariants \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O source.g.vcf.gz 

    # Then genotype both variant and invariant sites using the gvcf
    gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
        -R ${4} \
        -V source.g.vcf.gz \
        -L ${6} \
        -O calls.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals true \
        --only-output-calls-starting-in-intervals \
        --include-non-variant-sites true \
        --tmp-dir /tmp

     # Prepare annotation files for re-adding specific genotype fields to invariant sites

     # Subset gVCF to <NON_REF> sites, then convert to allsites vcf, and extract just chrom, pos, PL tags
     # NOTE any invariant sites with multiple alleles will be missed
    bcftools view -i 'N_ALT=1 && ALT="<NON_REF>"' source.g.vcf.gz \
        | bcftools convert --threads 4 --gvcf2vcf --fasta-ref  ${4} \
        | bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' \
        | sed 's/\t$//' \
        | bgzip > source_PL.tsv.gz
    tabix -s1 -b2 -e2 source_PL.tsv.gz
     
     # Subset gVCF to <NON_REF> sites, then convert to allsites vcf, and extract just chrom, pos, GQ tags
     # NOTE any invariant sites with multiple alleles will be missed
	bcftools view -i 'N_ALT=1 && ALT="<NON_REF>"' source.g.vcf.gz \
        | bcftools convert --threads 4 --gvcf2vcf --fasta-ref ${4} \
        | bcftools query -f '%CHROM\t%POS\t[%GQ\t]\n' \
        | sed 's/\t$//' \
        | bgzip > source_GQ.tsv.gz
    tabix -s1 -b2 -e2 source_GQ.tsv.gz

     # Subset called vcf to missing sites, then extract just chrom, pos, DP tags, then fabricate an AD column from the DP column 
    bcftools view -i 'ALT="."' calls.vcf.gz \
        | bcftools query -f '%CHROM\t%POS\t[%DP,0\t]\n' \
        | sed 's/\t$//' \
        | bgzip > source_AD.tsv.gz
    tabix -s1 -b2 -e2 source_AD.tsv.gz

     # Annotate the GQ, PL, AD tags for those <NON_REF> sites that lost it during variant calling
     # Replace missing alleles with <NON_REF> to ensure compatibility with CalculateGenotypePosteriors, otherwise it fails
     bcftools annotate -a source_GQ.tsv.gz -c CHROM,POS,FORMAT/GQ calls.vcf.gz -Ou \
        | bcftools annotate -a source_PL.tsv.gz -c CHROM,POS,FORMAT/PL -Ou \
        | bcftools annotate -a source_AD.tsv.gz -c CHROM,POS,FORMAT/AD -Ov \
        | awk 'BEGIN{OFS="\t"}
            /^#/ {print; next}
            { if($5==".") $5="<NON_REF>"; print }' \
        | bgzip > joint_called.vcf.gz
    tabix joint_called.vcf.gz
fi 

# Calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx${2}G" CalculateGenotypePosteriors \
    -V joint_called.vcf.gz \
    -L ${6} \
    -O joint_called_posterior.vcf.gz \
    --interval-merging-rule ALL \
    --merge-input-intervals true \
    --tmp-dir /tmp

# Convert <NON_REF> back to missing to ensure compatibility with filtering steps
    bcftools view joint_called_posterior.vcf.gz -Ov \
        | awk 'BEGIN{OFS="\t"}
            /^#/ {print; next}
            { if($5=="<NON_REF>") $5="."; print }' \
        | bgzip > ${5}.vcf.gz
    tabix ${5}.vcf.gz

# Remove temporary files
rm -f source.g.vcf.gz* calls.vcf.gz* *.tsv.gz* joint_called.vcf.gz* joint_called_posterior.vcf.gz*