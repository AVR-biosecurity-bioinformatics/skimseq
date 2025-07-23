#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = ref_genome

# Get prefix of vcf file for output name
prefix=$(echo ${2} | cut -d'.' -f1)

# Note - beagle output only works for polymorphic sites but indels and nonvariants are supported

# sample list
bcftools query -l ${2} > samples.txt

# Get Allele Depth (AD) matrix
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${2}  > ad_matrix.tsv

# build header
{ printf '#CHROM\tPOS'; while read s; do printf '\t%s' "$s"; done < samples.txt; printf '\n'; } > ph_calls_header.tsv

# draw pseudohaploid alleles
# Based on random sampling using allelic depth for each allele  
awk -F'\t' -v OFS='\t' -v seed=123 '
BEGIN{srand(seed)}
{
    chrom=$1;pos=$2;
    out=chrom OFS pos;
    for(i=5;i<=NF;i++){
        ad=$i;
        if(ad=="."||ad==""){
            call=".";
        }else{
            nAD=split(ad,counts,",");
            tot=0;
            for(j=1;j<=nAD;j++){c=counts[j]+0;if(c>0)tot+=c}
            if(tot==0){
                call=".";
            }else{
                r=rand()*tot;acc=0;
                for(j=1;j<=nAD;j++){acc+=counts[j];if(r<acc){call=j-1;break}}
            }
        }
        out=out OFS call;
    }
    print out;
}' ad_matrix.tsv > ph_calls_body.tsv

cat ph_calls_header.tsv ph_calls_body.tsv | bgzip > pseudohaploid_PH.tsv.gz
tabix -s1 -b2 -e2 pseudohaploid_PH.tsv.gz

# Add FORMAT header line
bcftools view -h ${2} | grep '^##' > hdr.txt
echo '##FORMAT=<ID=PH,Number=1,Type=String,Description="Pseudohaploid allele index sampled proportional to AD counts (0=REF,1=ALT1,...)">' >> hdr.txt
bcftools view -h ${2} | grep '^#CHROM' >> hdr.txt
bcftools reheader -h hdr.txt ${2} -o tmp.reheader.vcf.gz
bcftools index -t tmp.reheader.vcf.gz

# Annotate PH tag onto existing VCF
bcftools annotate -a pseudohaploid_PH.tsv.gz -c CHROM,POS,FORMAT/PH -Oz -o raw.withPH.vcf.gz tmp.reheader.vcf.gz
bcftools index -t raw.withPH.vcf.gz

# Replace GT tag with new PH tag
# 1 set 0/0 where PH==0
# 2 set 1/1 where PH==1
# 3 set ./., where PH is missing
bcftools +setGT raw.withPH.vcf.gz -Ou -- -t q -n c:0/0 -i 'FMT/PH=="0"' \
    | bcftools +setGT -Ou -- -t q -n c:1/1 -i 'FMT/PH=="1"' \
    | bcftools +setGT -Oz -o ${prefix}_pseudohaploid.vcf.gz -- -t q -n . -i 'FMT/PH=="."' 

bcftools index -t ${prefix}_pseudohaploid.vcf.gz

# Output genotype matrix (angsd style)

# Create header
{
  printf "chr\tpos\tmajor\tminor";
  bcftools query -l ${prefix}_pseudohaploid.vcf.gz | while read s; do printf '\t%s' "$s"; done
  printf '\n'
} > ${prefix}_pseudohaploid.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PH]\n' ${prefix}_pseudohaploid.vcf.gz \
    | awk -F'\t' 'BEGIN{OFS="\t"}
    {
    chr=$1; pos=$2; ref=$3; alt=$4;
    printf "%s\t%s\t%s\t%s", chr, pos, ref, alt

    for (i=5; i<=NF; i++) {
        if ($i=="0")      out=0;       # REF
        else if ($i=="1") out=1;       # ALT
        else              out=-1;      # missing or unexpected
        printf "\t%s", out
    }
    printf "\n"
    }' >> ${prefix}_pseudohaploid.tsv