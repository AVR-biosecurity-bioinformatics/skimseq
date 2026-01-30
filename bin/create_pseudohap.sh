#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = ref_genome

VCF=${2}

# sample list
bcftools query -l ${VCF} > samples.txt

# build header
{ printf '#CHROM\tPOS'; while read s; do printf '\t%s' "$s"; done < samples.txt; printf '\n'; } > ph_VCF_header.tsv

#  Get Allele Depth (AD) matrix, and draw pseudohaploid alleles based on random sampling using allelic depth for each allele  
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${VCF}  | awk -F'\t' -v OFS='\t' -v seed=123 '
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
                r=int(rand()*tot);acc=0;
                for(j=1;j<=nAD;j++){acc+=counts[j];if(r<acc){call=j-1;break}}
            }
        }
        out=out OFS call;
    }
    print out;
}' > ph_VCF_body.tsv

cat ph_VCF_header.tsv ph_VCF_body.tsv | bgzip > pseudohaploid_PH.tsv.gz
tabix -s1 -b2 -e2 pseudohaploid_PH.tsv.gz

# Add FORMAT header line
bcftools view -h ${VCF} | grep '^##' > hdr.txt
echo '##FORMAT=<ID=PH,Number=1,Type=String,Description="Pseudohaploid allele index sampled proportional to AD counts (0=REF,1=ALT1,...)">' >> hdr.txt
bcftools view -h ${VCF} | grep '^#CHROM' >> hdr.txt
bcftools reheader -h hdr.txt ${VCF} -o tmp.reheader.vcf.gz
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
    | bcftools +setGT -Oz -o pseudohaploid.vcf.gz -- -t q -n . -i 'FMT/PH=="."' 

bcftools index -t $pseudohaploid.vcf.gz