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

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
GENOMICSDB="${3}"
REF="${4}"
IHASH="${5}"
INTERVAL_BED="${6}"
EXCLUDE_BED="${7}"

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${MEM_GB} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
# Send stderr to log file for profiling
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}G" GenotypeGVCFs \
    -R "${REF}" \
    -V gendb://${GENOMICSDB} \
    -L "${INTERVAL_BED}" \
    -O genotyped.vcf.gz \
    --exclude-intervals "${EXCLUDE_BED}" \
    --interval-exclusion-padding "${EXCLUDE_PAD}" \
    --include-non-variant-sites "${OUTPUT_INVARIANT}" \
    --interval-merging-rule ALL \
    --merge-input-intervals \
    --variant-output-filtering STARTS_IN \
    --max-alternate-alleles "${MAX_ALTERNATE}" \
    --genomicsdb-max-alternate-alleles "${GENOMICSDB_MAX_ALTERNATE}" \
    -ploidy "${PLOIDY}" \
    --heterozygosity "${HET}" \
    --heterozygosity-stdev "${HET_SD}" \
    --indel-heterozygosity "${INDEL_HET}" \
    --tmp-dir /tmp \
    --genomicsdb-shared-posixfs-optimizations true \
    2> >(tee -a ${IHASH}.stderr.log >&2)

# Add additional annotations to interval VCF file - To be used for filtering

# Add minor alelle count (MAC) info tag
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
bcftools annotate --threads ${CPUS} -h MAC.hdr -a MAC.tsv.gz -c CHROM,POS,INFO/MAC -Ou genotyped.vcf.gz \
    | bcftools +setGT -- -t q -n . -i 'FMT/DP=0' \
    | bcftools +fill-tags -- -t AC_Hom,AC_Het,AC_Hemi,MAF,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \
    | bcftools +tag2tag -- --PL-to-GL \
    | bcftools annotate --threads ${CPUS} --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Ou \
    | bcftools sort -Oz9 -o ${IHASH}.vcf.gz 

# Reindex outpu
bcftools index -t ${IHASH}.vcf.gz

# Clean up
rm genotyped.vcf.gz* *.hdr *.tsv.gz