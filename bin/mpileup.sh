#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref_genome
# $4 = interval hash
# $5 = interval_bed

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
IHASH="${4}"
INTERVAL_BED="${5}"

# -----------------------------
# Pre-filter reads using samtools
# -----------------------------

# Set up filtering expressions for samtools

# Min aligned bases (excluding softclips)
ALN_CLAUSE="(qlen - sclen) >= ${MIN_ALIGNED_LENGTH}"

# Min and max fragment (insert size) length
TLEN_CLAUSE="( (tlen >= ${MIN_FRAGMENT_LENGTH} && tlen <= ${MAX_FRAGMENT_LENGTH}) || (tlen <= -${MIN_FRAGMENT_LENGTH} && tlen >= -${MAX_FRAGMENT_LENGTH}) )"

# Read flags filter (make DUP conditional)
if [[ "${RMDUP}" == "false" ]]; then
  FLAGS_CLAUSE='(!flag.unmap && !flag.secondary && !flag.supplementary  && flag.proper_pair )'
else
  FLAGS_CLAUSE='(!flag.unmap && !flag.secondary && !flag.supplementary  && flag.proper_pair && !flag.dup)'
fi

# Create joint filtering expression for samtools
EXPR="${FLAGS_CLAUSE} && ${ALN_CLAUSE} && ${TLEN_CLAUSE}"

export REF INTERVAL_BED EXPR

# Choose jobs and cpus for GNU parallel
THREADS_PER_JOB=2
JOBS=$(( CPUS / THREADS_PER_JOB )) # how many CRAMs in parallel
(( JOBS < 1 )) && JOBS=1

# Subset and filter CRAMs in parallel
parallel --jobs "${JOBS}" --line-buffer '
  cram={}
  out="${cram%.cram}.filt.cram"

  samtools view -@ '"${THREADS_PER_JOB}"' -T "${REF}" \
    --regions-file "${INTERVAL_BED}" \
    -e "${EXPR}" \
    -O cram -o "${out}" "${cram}"

  samtools index -@ '"${THREADS_PER_JOB}"' "${out}"
  echo "${out}"
' :::: cram.list > cram.filtered.list

# -----------------------------
# Variant calling using pre-filtered crams
# -----------------------------

# set up bcftools filter flags
if [[ "${RMDUP}" == "false" ]]; then
  FILTER_FLAGS="--ns DUP -G UNMAP,SECONDARY,QCFAIL"
else
  FILTER_FLAGS="--ns UNMAP,SECONDARY,QCFAIL,DUP"
fi

if [[ "${OUTPUT_INVARIANT}" == "false" ]]; then
  VARIANTS_ONLY="--variants-only"
else
  VARIANTS_ONLY=""
fi

# Call variants with mpleup
bcftools mpileup \
    --threads ${CPUS} \
    --bam-list cram.filtered.list \
    --max-depth ${MAXDEPTH} \
    --fasta-ref ${REF} \
    --min-BQ ${MINBQ} \
    --min-MQ ${MINMQ} \
    --regions-file ${INTERVAL_BED} \
    ${FILTER_FLAGS} \
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
    --prior ${MUTATION_RATE}

# -----------------------------
# Add additional annotations to be used for filtering
# -----------------------------

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
    | bcftools +fill-tags -- -t MAF,ExcHet,HWE,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \
    | bcftools annotate --threads ${CPUS} --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -Ou \
    | bcftools sort -Oz9 -o ${IHASH}.vcf.gz 

#| bcftools +tag2tag -- --PL-to-GL \

# Reindex output
bcftools index -t ${IHASH}.vcf.gz

# Clean up temporary files
xargs -r -d '\n' rm -f < <(awk '{print $0; print $0 ".crai"}' cram.filtered.list)
rm -f cram.filtered.list
rm genotyped.vcf.gz* *.hdr *.tsv.gz