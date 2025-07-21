#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem

# get list of .vcf files in directory
ls *.vcf.gz > vcf.list

# merge all .vcfs together
gatk --java-options "-Xmx${2}G" MergeVcfs \
    -I vcf.list \
    -O merged_tmp.vcf.gz

# Add additional info tags to be used for filtering
bcftools +fill-tags --threads ${1} merged_tmp.vcf.gz -o merged_tmp2.vcf.gz -- -t AC,AN,MAF,F_MISSING,NS

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' merged_tmp2.vcf.gz \
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
     -O z -o merged.vcf.gz \
     merged_tmp2.vcf.gz

# reindex the output file
gatk --java-options "-Xmx${2}G" IndexFeatureFile \
    -I merged.vcf.gz 

# Remove the temporary files
rm -f merged_tmp.vcf.gz* merged_tmp2.vcf.gz* MAC.tsv.gz*