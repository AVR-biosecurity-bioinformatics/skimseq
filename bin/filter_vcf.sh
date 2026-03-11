#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed
# $6 = interval_hash
# $7 = sample_groups_tsv
# $8 = missing data summary

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# TODO: Break out vcf masks into individual components (i.e. Genmap, longdust, etc)

# Map variant type to bcftools selectors
# TODO: select variants on output of joint genotype
case "${4}" in
  snp)       TYPE_ARGS="-v snps   -m2 -M2 -e 'ALT=\"*\"'";;   # biallelic SNPs, drop star alleles
  indel)     TYPE_ARGS="-v indels -m2 -M2";;                  # biallelic INDELs
  invariant) TYPE_ARGS="-v ref";;                              # reference-only sites (if present)
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

# Find samples above the missing fraction filter
awk -v thr="$SAMPLE_MAX_MISSING" 'NR==1 {next} $4!="NA" && ($4+0) < thr {print $1}' "${8}" > samples_to_keep.txt

# Create sample_groups.tsv:
# first keep only samples in samples_to_keep.txt
# then drop populations with fewer than MIN_SAMPLES_PER_POP retained samples
awk -v n="$MIN_SAMPLES_PER_POP" '
    BEGIN { FS=OFS="\t" }
    NR==FNR {
        keep[$1] = 1
        next
    }
    ($1 in keep) {
        count[$2]++
        sample[NR] = $1
        pop[NR]    = $2
    }
    END {
        for (i = 1; i <= NR; i++) {
            if (sample[i] != "" && count[pop[i]] >= n) {
                print sample[i], pop[i]
            }
        }
    }
' samples_to_keep.txt "${7}" > sample_groups.tsv

#TODO: need to rename any samples with 2 letter names

# number of pops (one vcf per pop)
N_POPS=$(cut -f2 sample_groups.tsv | tr ',' '\n' | sort -u | sed '/^$/d' | wc -l)

# percentage-based threshold: Tp = ceil(PERC * n_pops / 100)
# (use ceil so 1% of 3 pops => 1, not 0; change to floor if you prefer)
PERC_N=$(awk -v p="$PERC_POPS_FAILING" -v n="$N_POPS" '
  BEGIN{
    if(p<=0 || n<=0){print 0; exit}
    x = p*n
    t = int(x)
    if(t < x) t++
    print t
  }'
)

# Need defaults here because the exporting of kv pairs unsets anything thats -1
N_POPS_FAILING_VAL="${N_POPS_FAILING:--1}"
PERC_N_VAL="${PERC_N:--1}"

# choose effective minimum number of populations
if (( N_POPS_FAILING_VAL < 0 && PERC_N_VAL < 0 )); then
  MIN_POPS=-1
elif (( N_POPS_FAILING_VAL < 0 )); then
  MIN_POPS=$PERC_N_VAL
elif (( PERC_N_VAL < 0 )); then
  MIN_POPS=$N_POPS_FAILING_VAL
elif (( PERC_N_VAL > N_POPS_FAILING_VAL )); then
  MIN_POPS=$PERC_N_VAL
else
  MIN_POPS=$N_POPS_FAILING_VAL
fi

# Helper function to make per-population filter expressions
make_pop_count_expr() {
  local tag="$1"      # e.g. MAF
  local op="$2"        # >, >=, <, <=
  local thr="$3"      # e.g. 0.0005
  local min_n="$4"    # e.g. 2
  local pops_file="$5"

  cut -f2 "$pops_file" \
    | tr ',' '\n' \
    | sort -u \
    | awk -v tag="$tag" -v op="$op" -v thr="$thr" -v n="$min_n" '
        BEGIN { first=1 }
        {
          if (!first) printf " + "
          printf "(INFO/%s_%s%s%s)", tag, $1, op, thr
          first=0
        }
        END { printf " < %s", n }
      '
}

# Per-pop fail expressions:
# fail if fewer than MIN_POPS populations have tag meeting the criterion
pop_eh_expr=$(make_pop_count_expr "ExcHet" ">=" "${EH:--1}" "$MIN_POPS" sample_groups.tsv)
pop_hwe_expr=$(make_pop_count_expr "HWE" ">=" "${HWE:--1}" "$MIN_POPS" sample_groups.tsv)
pop_maf_expr=$(make_pop_count_expr "MAF" ">=" "${MAF:-0}" "$MIN_POPS" sample_groups.tsv)
pop_mac_expr=$(make_pop_count_expr "MAC" ">=" "${MAC:-0}" "$MIN_POPS" sample_groups.tsv)
pop_ns_expr=$(make_pop_count_expr "NS" ">=" "${NS:-0}" "$MIN_POPS" sample_groups.tsv)
pop_cr_expr=$(make_pop_count_expr "CR" ">=" "${CR:-0}" "$MIN_POPS" sample_groups.tsv)

# Subset to target variant class and just samples above missing data filter
# Then add global annotations using fill-tags
# Then add per-pop annotations  using fill-tags
# Note MAC is calculated from MAF (7 decimal precision), this could cause rounding for very large cohorts (i.e. 100k+)
bcftools view --threads ${1} ${TYPE_ARGS} -S samples_to_keep.txt -Ou "${3}" \
  | bcftools +setGT -Ou -- \
    -t q \
    -n . \
    -i "FORMAT/GQ < ${GQ:-0} || FORMAT/DP < ${gtDPmin:-0} || FORMAT/DP > ${gtDPmax:-999999999}" \
  | bcftools +fill-tags -Ou - -- \
    -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING,MAC=int(MAF*AN)' \
  | bcftools +fill-tags -Ou - -- \
    -S sample_groups.tsv \
    -t 'AC,AN,NS,MAF,HWE,ExcHet,CR:1=1-F_MISSING,MAC=int(MAF*AN)' \
  | bcftools filter -Ou -s MASK_FAIL       -m+ -M vcf_masks.bed \
  | bcftools filter -Ou -s QUAL_FAIL       -m+ -e "QUAL < ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s DP_FAIL         -m+ -e "INFO/DP < ${DPmin:-0} || INFO/DP < ${DPlower:-0} || INFO/DP > ${DPupper:-999999999}" \
  | bcftools filter -Ou -s DIST_INDEL_FAIL -m+ -e "INFO/DIST_INDEL < ${DIST_INDEL:--999999999}" \
  | bcftools filter -Ou -s EH_FAIL         -m+ -e "$pop_eh_expr" \
  | bcftools filter -Ou -s HWE_FAIL        -m+ -e "$pop_hwe_expr" \
  | bcftools filter -Ou -s MAF_FAIL        -m+ -e "$pop_maf_expr" \
  | bcftools filter -Ou -s MAC_FAIL        -m+ -e "$pop_mac_expr" \
  | bcftools filter -Ou -s NS_FAIL         -m+ -e "$pop_ns_expr" \
  | bcftools filter -Ou -s CR_FAIL         -m+ -e "$pop_cr_expr" \
  | bcftools view --threads ${1} -Ob -o tmp.bcf

# Drop failing genotypes to create filtered vcf file (main output)
bcftools view --threads "${1}" -f PASS -Oz9 -o ${4}.${6}.filt.vcf.gz tmp.bcf
bcftools index --threads ${1} -t ${4}.${6}.filt.vcf.gz

# Create filtered sitelist file by dropping non-passing variants and genotypes
bcftools view --threads "${1}" -G -f PASS -Oz9 -o ${4}.${6}.sitelist.vcf.gz tmp.bcf
bcftools index --threads ${1} -t ${4}.${6}.sitelist.vcf.gz

# Create sites only tagged file for QC histograms
bcftools view --threads "${1}" -G -Oz9 -o ${4}.${6}.tagged.vcf.gz tmp.bcf
bcftools index --threads ${1} -t ${4}.${6}.tagged.vcf.gz

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n ${4}.${6}.filt.vcf.gz | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${6}.counts"