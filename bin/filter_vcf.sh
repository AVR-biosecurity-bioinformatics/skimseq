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
awk -v thr="$SAMPLE_MAX_MISSING" 'NR==1 {next} $4!="NA" && ($4+0) < thr {print $1}' "${8}" > ${4}.${6}.samples.txt

# Create sample_groups.tsv:
# first keep only samples in ${4}.${6}.samples.txt
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
' ${4}.${6}.samples.txt "${7}" > sample_groups.tsv

# number of pops (one vcf per pop)
N_POPS=$(cut -f2 sample_groups.tsv | tr ',' '\n' | sort -u | sed '/^$/d' | wc -l)

# Convert the proportional population threshold into an absolute count using ceil
# (e.g. 0.9 * 5 pops = 4.5 -> 5 pops must pass)
PERC_N=$(awk -v p="${PERC_POPS_FAILING:-}" -v n="$N_POPS" '
  BEGIN{
    if (p == "" || p == "NA" || p == "na" || n <= 0) { print ""; exit }
    x = p * n
    t = int(x)
    if (t < x) t++
    print t
  }'
)

# Use whichever population-pass threshold is more stringent:
# - N_POPS_FAILING converted externally to a minimum pass requirement
# - PERC_N derived from the proportion threshold above
if [[ -z "${N_POPS_FAILING:-}" && -z "$PERC_N" ]]; then
  MIN_POPS=""
elif [[ -z "${N_POPS_FAILING:-}" ]]; then
  MIN_POPS="$PERC_N"
elif [[ -z "$PERC_N" ]]; then
  MIN_POPS="$N_POPS_FAILING"
elif (( PERC_N > N_POPS_FAILING )); then
  MIN_POPS="$PERC_N"
else
  MIN_POPS="$N_POPS_FAILING"
fi

# Fail when greater or equal to MIN_POPS
MIN_PASS=$(( N_POPS - MIN_POPS + 1 ))
echo "MIN_PASS=$MIN_PASS"

# Helper functions for per-population fail tags
# make_pass_tags() determines whether the annotated value for a pop passes the per-pop filter threshold
# This generates new per-pop tags like: MAF_PASS_Pop1:1=int(MAF_Pop1<0.05)
make_pass_tags() {
  local base="$1"      # e.g. MAF, NS, HWE, ExcHet
  local op="$2"        # pass operator, usually >=
  local thr="$3"       # threshold
  local pops_file="$4"

  # disabled
  if [[ -z "$thr" ]]; then
    echo ""
    return
  fi

  cut -f2 "$pops_file" \
    | tr ',' '\n' \
    | sed '/^$/d' \
    | sort -u \
    | awk -v base="$base" -v op="$op" -v thr="$thr" '
        BEGIN { first=1 }
        {
          if (!first) printf ","
          printf "%s_PASS_%s:1=int(%s_%s%s%s)", base, $1, base, $1, op, thr
          first=0
        }'
}


# make_npass_tag() counts the number of passing populations and adds a new annotation
# This generates a single count tag like: NFAIL_MAF:1=int(MAF_FAIL_Pop1+MAF_FAIL_Pop2+...)
make_npass_tag() {
  local base="$1"
  local pops_file="$2"

  local sum_expr
  sum_expr=$(
    cut -f2 "$pops_file" \
      | tr ',' '\n' \
      | sed '/^$/d' \
      | sort -u \
      | awk -v base="$base" '
          BEGIN { first=1 }
          {
            if (!first) printf "+"
            printf "%s_PASS_%s", base, $1
            first=0
          }'
  )

  printf 'NPASS_%s:1=int(%s)\n' "$base" "$sum_expr"
}

# Join non-empty strings with commas
join_tags() {
  local out=""
  local x
  for x in "$@"; do
    [[ -z "$x" ]] && continue
    [[ -n "$out" ]] && out+=","
    out+="$x"
  done
  printf '%s\n' "$out"
}

# Set up pass flags for each filter and pop and concatenate together into single string for 1 pass annotation
POP_PASS_TAGS=$(join_tags \
"$(make_pass_tags "ExcHet" ">=" "${EH_POP:-}"  sample_groups.tsv)" \
"$(make_pass_tags "HWE"    ">=" "${HWE_POP:-}" sample_groups.tsv)" \
"$(make_pass_tags "MAF"    ">=" "${MAF_POP:-}" sample_groups.tsv)" \
"$(make_pass_tags "NS"      ">=" "${NS_POP:-}"  sample_groups.tsv)"
)

# Set up pass counting tags for each filter and concatenate together into single string for 1 pass annotation
NPASS_TAGS=$(join_tags \
  "$(make_npass_tag "ExcHet" sample_groups.tsv)" \
  "$(make_npass_tag "HWE"    sample_groups.tsv)" \
  "$(make_npass_tag "MAF"    sample_groups.tsv)" \
  "$(make_npass_tag "NS"      sample_groups.tsv)"
)

echo "POP_PASS_TAGS:"
echo $POP_PASS_TAGS

echo "NPASS_TAGS:"
echo $NPASS_TAGS

# Subset to target variant class and just samples above missing data filter
# Then add global annotations using fill-tags
# Then add per-pop threshold annotations using fill-tags
# Then add per-population pass flags, using the per-population threshold annotations using fill-tags
# Then count the number of populations that pass
# Note MAC is calculated from MAF (7 decimal precision), this could cause rounding for very large cohorts (i.e. 100k+)
# NOTE: bcftools +fill-tags breaks with any samples that have 2 letter names
bcftools view --threads ${1} ${TYPE_ARGS} -S ${4}.${6}.samples.txt -Ou "${3}" \
  | bcftools +setGT -Ou -- \
    -t q \
    -n . \
    -i "FORMAT/GQ < ${GQ:-0} | FORMAT/DP < ${gtDPmin:-0} | FORMAT/DP > ${gtDPmax:-999999999}" \
  | bcftools +fill-tags -Ou - -- \
    -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING' \
  | bcftools +fill-tags -Ou - -- \
    -S sample_groups.tsv \
    -t 'NS,MAF,HWE,ExcHet' \
  | bcftools +fill-tags -Ou - -- \
    -t "$POP_PASS_TAGS" \
  | bcftools +fill-tags -Ou - -- \
    -t "$NPASS_TAGS" \
  | bcftools filter -Ou -s MASK_FAIL       -m+ -M vcf_masks.bed \
  | bcftools filter -Ou -s QUAL_FAIL       -m+ -e "QUAL < ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s DP_FAIL         -m+ -e "INFO/DP < ${DPmin:-0} || INFO/DP < ${DPlower:-0} || INFO/DP > ${DPupper:-999999999}" \
  | bcftools filter -Ou -s DIST_INDEL_FAIL -m+ -e "INFO/DIST_INDEL < ${DIST_INDEL:--999999999}" \
  | bcftools filter -Ou -s EH_FAIL         -m+ -e "INFO/ExcHet < ${EH_GLOBAL:--1}" \
  | bcftools filter -Ou -s HWE_FAIL        -m+ -e "INFO/HWE < ${HWE_GLOBAL:--1}" \
  | bcftools filter -Ou -s MAF_FAIL        -m+ -e "INFO/MAF < ${MAF_GLOBAL:-0}" \
  | bcftools filter -Ou -s NS_FAIL         -m+ -e "INFO/NS < ${NS_GLOBAL:-0}" \
  | bcftools filter -Ou -s CR_FAIL         -m+ -e "INFO/CR < ${CR_GLOBAL:-0}" \
  | bcftools filter -Ou -s POP_EH_FAIL     -m+ -e "${MIN_PASS:+INFO/NPASS_ExcHet < ${MIN_PASS}}" \
  | bcftools filter -Ou -s POP_HWE_FAIL    -m+ -e "${MIN_PASS:+INFO/NPASS_HWE < ${MIN_PASS}}" \
  | bcftools filter -Ou -s POP_MAF_FAIL    -m+ -e "${MIN_PASS:+INFO/NPASS_MAF < ${MIN_PASS}}" \
  | bcftools filter -Ou -s POP_NS_FAIL     -m+ -e "${MIN_PASS:+INFO/NPASS_NS < ${MIN_PASS}}" \
  | bcftools view --threads ${1} -Ob -o tmp.bcf

# Drop failing sites to create filtered vcf file (main output)
bcftools view --threads "${1}" -f PASS -Ou tmp.bcf \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' -Oz9 -o ${4}.${6}.filt.vcf.gz 
bcftools index --threads ${1} -t ${4}.${6}.filt.vcf.gz

# Create filtered sitelist file by dropping non-passing variants and genotypes
bcftools view --threads "${1}" -G -f PASS -Ou tmp.bcf \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' -Oz9 -o ${4}.${6}.sitelist.vcf.gz 
bcftools index --threads ${1} -t ${4}.${6}.sitelist.vcf.gz

# Create sites only tagged file for QC histograms
bcftools view --threads "${1}" -G -Oz9 -o ${4}.${6}.tagged.vcf.gz tmp.bcf
bcftools index --threads ${1} -t ${4}.${6}.tagged.vcf.gz

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n ${4}.${6}.filt.vcf.gz | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.${6}.counts"