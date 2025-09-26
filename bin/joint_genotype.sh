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

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
if [[ "${9}" == "false" ]]; then
    # joint genotype variant sites only
    gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" GenotypeGVCFs \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O joint_called.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals \
        --only-output-calls-starting-in-intervals \
        --max-alternate-alleles 6 \
        --genomicsdb-max-alternate-alleles 10 \
        --tmp-dir /tmp \
        --genomicsdb-shared-posixfs-optimizations true

elif [[ "${9}" == "true" ]]; then
    # Joint genotype both variant and invariant
    # This requires some custom code to re-add missing genotpye fields for compatibility with later steps

    # genotype both variant and invariant sites 
    gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g"  GenotypeGVCFs \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O calls.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals true \
        --only-output-calls-starting-in-intervals \
        --max-alternate-alleles 6 \
        --genomicsdb-max-alternate-alleles 10 \
        --tmp-dir /tmp \
        --genomicsdb-shared-posixfs-optimizations true

    # Get the sites as a GVCF as well to transfer the annotations over
    gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g"  SelectVariants \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O source.g.vcf.gz 

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
else 
    echo "output_invariant must be true or false"
    exit 1
fi 

# Calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" CalculateGenotypePosteriors \
    -V joint_called.vcf.gz \
    -L ${6} \
    -O joint_called_posterior.vcf.gz \
    --interval-merging-rule ALL \
    --tmp-dir /tmp

# Convert any <NON_REF> back to missing to ensure compatibility with filtering steps
bcftools view joint_called_posterior.vcf.gz -Ov \
    | awk 'BEGIN{OFS="\t"}
        /^#/ {print; next}
        { if($5=="<NON_REF>") $5="."; print }' \
    | bgzip > tmp.vcf.gz
tabix tmp.vcf.gz

# Add additional info tags to be used for filtering
bcftools +fill-tags --threads ${1} tmp.vcf.gz -o tmp2.vcf.gz -- -t AC,AN,MAF,F_MISSING,NS

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp2.vcf.gz \
| awk 'BEGIN{OFS="\t"}
       {
         split($3,ac,",")          # AC is comma‑separated if multi‑allelic
         mac=$4                    # start with AN
         refCount = $4             # will be AN - sum(AC)
         for(i in ac){refCount-=ac[i]; mac=(ac[i]<mac?ac[i]:mac)}
         mac=(refCount<mac?refCount:mac)
         print $1,$2,mac
       }'                         \
| bgzip > MAC.tsv.gz
tabix -s1 -b2 -e2 MAC.tsv.gz

# Create VCF header line for MAC filter
cat > add_MAC.hdr <<'EOF'
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count (minimum of each ALT AC and reference allele count)">
EOF

# Annotate the vcf with INFO/MAC
bcftools annotate \
     -h add_MAC.hdr \
     -a MAC.tsv.gz \
     -c CHROM,POS,INFO/MAC \
     -O z -o ${5}.vcf.gz \
     tmp2.vcf.gz

# reindex the output file
tabix ${5}.vcf.gz

# Remove temporary files
rm -f source.g.vcf.gz* calls.vcf.gz* *.tsv.gz* joint_called.vcf.gz* joint_called_posterior.vcf.gz* tmp.vcf.gz* tmp2.vcf.gz*