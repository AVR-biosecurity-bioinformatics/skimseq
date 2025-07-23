#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf

# Add additional info tags to be used for filtering
bcftools +fill-tags --threads ${1} ${3} -o tmp.vcf.gz -- -t AC,AN,MAF,F_MISSING,NS

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp.vcf.gz \
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
     -O z -o annot.vcf.gz \
     tmp.vcf.gz

# reindex the output file
gatk --java-options "-Xmx${2}G" IndexFeatureFile \
    -I annot.vcf.gz

# Remove the temporary files
rm -f tmp.vcf.gz* *.tsv.gz*