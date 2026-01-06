#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref_genome
# $4 = interval hash
# $5 = interval_bed
# $6 = exclude_bed
# $7 = cram

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
IHASH="${4}"
INTERVAL_BED="${5}"
EXCLUDE_BED="${6}"

# Create new bed file, subtracting exclude_bed+padding

# parse filtering options as flags
if [[ ${RMDUP} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    NS="--ns DUP -G UNMAP,SECONDARY,QCFAIL"
else 
    NS="--ns UNMAP,SECONDARY,QCFAIL,DUP"
fi

if [[ ${OUTPUT_INVARIANT} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    VARIANTS_ONLY="--variants-only"
else 
    VARIANTS_ONLY=""
fi

# Create list of crams to be processed
ls *.cram | grep -v '.crai' | sort > cram.list

# TODO - need to handle whether duplicates should be removed
# TODO - need to toggle variants only
bcftools mpileup \
    --threads ${CPUS} \
    --bam-list cram.list \
    --max-depth 250 \
    --fasta-ref ${REF} \
    --min-BQ ${MINBQ} \
    --min-MQ ${MINMQ} \
    --regions-file ${INTERVAL_BED} \
    ${NS} \
    --annotate FORMAT/DP,INFO/AD,FORMAT/DP,FORMAT/SP \
    --indels-cns \
    --indel-size 110 \
    | bcftools call \
    -Ob \
    -o genotyped.vcf.gz \
    -a FORMAT/GP,FORMAT/GQ,INFO/PV4 \
    --ploidy ${PLOIDY} \
    ${VARIANTS_ONLY} \
    --multiallelic-caller \
    --prior 1.1e-3 \
    --pval-threshold 0.5
   

# Add additional annotations to interval VCF file - To be used for filtering

# Create annotation table for minor allele count (MAC)
# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' genotyped.vcf.gz \
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
cat > MAC.hdr <<'EOF'
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count (minimum of each ALT AC and reference allele count)">
EOF

# Create annotation table for distance to closest indel
bcftools view -v indels -Ou genotyped.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%REF\n' \
| awk 'BEGIN{OFS="\t"}
       { s=$2-1; e=s+length($3); if(e<=s) e=s+1; print $1,s,e }' \
| LC_ALL=C sort -k1,1 -k2,2n > indels.bed

# Get a list of all sites
bcftools query -f '%CHROM\t%POS\n' genotyped.vcf.gz \
| awk 'BEGIN{OFS="\t"}{ print $1, $2-1, $2, $2 }' \
| LC_ALL=C sort -k1,1 -k2,2n > sites.bed

# Find distance to closest indel. No indel returns -1
bedtools closest -a sites.bed -b indels.bed -d -t first \
| awk 'BEGIN{OFS="\t"} {print $1,$4,$NF}' \
| bgzip > dist_to_indel.tsv.gz

tabix -s1 -b2 -e2 dist_to_indel.tsv.gz

cat > dist_to_indel.hdr <<'EOF'
##INFO=<ID=DIST_INDEL,Number=1,Type=Integer,Description="Distance in bp to closest indel (0 if overlaps indel reference span)">
EOF

# Add annotations to vcf
bcftools annotate --threads ${CPUS} -h MAC.hdr -a MAC.tsv.gz -c CHROM,POS,INFO/MAC -Ou genotyped.vcf.gz \
    | bcftools annotate --threads ${CPUS} -h dist_to_indel.hdr -a dist_to_indel.tsv.gz -c CHROM,POS,INFO/DIST_INDEL -Ou \
    | bcftools +setGT -- -t q -n . -i 'FMT/DP=0' \
    | bcftools +fill-tags -- -t AC_Hom,AC_Het,AC_Hemi,MAF,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \
    | bcftools +tag2tag -- --PL-to-GL \
    | bcftools annotate --threads ${CPUS} --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -Ou \
    | bcftools sort -Oz9 -o ${IHASH}.vcf.gz 

# Reindex outpu
bcftools index -t ${IHASH}.vcf.gz

# Clean up
rm genotyped.vcf.gz* *.hdr *.tsv.gz