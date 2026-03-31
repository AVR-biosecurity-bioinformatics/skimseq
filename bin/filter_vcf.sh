#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = interval_hash
# $5 = mask_bed
# $6 = sample_groups_tsv
# $7 = missing data summary

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# TODO: Break out vcf masks into individual components (i.e. Genmap, longdust, etc)

# Find samples above the missing fraction filter
awk -v thr="$SAMPLE_MAX_MISSING" 'NR==1 {next} $4!="NA" && ($4+0) < thr {print $1}' "${7}" > ${4}.samples.txt

# Create sample_groups.tsv:
# first keep only samples in ${4}.samples.txt
# then drop populations with fewer than MIN_SAMPLES_PER_POP retained samples
awk -v n="$POPULATION_MIN_SAMPLES_PER_POP" '
    BEGIN { FS=OFS="\t" }
    NR==FNR {
        keep[$1] = 1
        next
    }
    ($1 in keep) {
        pop_name = $2
        gsub(/[[:space:]_]+/, "", pop_name)

        count[pop_name]++
        sample[NR] = $1
        pop[NR]    = pop_name
    }
    END {
        for (i = 1; i <= NR; i++) {
            if (sample[i] != "" && count[pop[i]] >= n) {
                print sample[i], pop[i]
            }
        }
    }
' "${4}.samples.txt" "${6}" > sample_groups.tsv

# number of pops (one vcf per pop)
N_POPS=$(cut -f2 sample_groups.tsv | tr ',' '\n' | sort -u | sed '/^$/d' | wc -l)

# Convert the proportional population threshold into an absolute count using ceil
# (e.g. 0.9 * 5 pops = 4.5 -> 5 pops must pass)
PERC_N=$(awk -v p="${POPULATION_PERC_POPS_FAILING:-}" -v n="$N_POPS" '
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
if [[ -z "${POPULATION_N_POPS_FAILING:-}" && -z "$PERC_N" ]]; then
  MIN_POPS=""
elif [[ -z "${POPULATION_N_POPS_FAILING:-}" ]]; then
  MIN_POPS="$PERC_N"
elif [[ -z "$PERC_N" ]]; then
  MIN_POPS="$POPULATION_N_POPS_FAILING"
elif (( PERC_N > POPULATION_N_POPS_FAILING )); then
  MIN_POPS="$PERC_N"
else
  MIN_POPS="$POPULATION_N_POPS_FAILING"
fi

# Fail when greater or equal to MIN_POPS
MIN_PASS=""
if [[ -n "${MIN_POPS:-}" ]]; then
  MIN_PASS=$(( N_POPS - MIN_POPS + 1 ))
fi

# Helper functions for per-population fail tags
# make_typed_pass_tags() determines whether the annotated value for a pop passes the per-pop filter threshold, in a variant type specfiic way
# This generates new per-pop tags like: MAF_PASS_Pop1:1=int(MAF_Pop1<0.05)
make_typed_pass_tags() {
  local base="$1"      # e.g. MAF
  local op="$2"        # usually >=
  local thr="$3"       # threshold
  local suffix="$4"    # e.g. SNP / INDEL / INVARIANT
  local pops_file="$5"

  cut -f2 "$pops_file" \
    | tr ',' '\n' \
    | sed '/^$/d' \
    | sort -u \
    | awk -v base="$base" -v op="$op" -v thr="$thr" -v suffix="$suffix" '
        BEGIN { first=1 }
        {
          if (!first) printf ","

          if (thr == "" || thr == "NA" || thr == "na") {
            printf "%s_PASS_%s_%s:1=1", base, suffix, $1
          } else {
            printf "%s_PASS_%s_%s:1=int(%s_%s%s%s)", base, suffix, $1, base, $1, op, thr
          }

          first=0
        }'
}

# make_typed_npass_tag() counts the number of passing populations and adds a new annotation, in a variant type specific way
# This generates a single count tag like: NFAIL_MAF:1=int(MAF_FAIL_Pop1+MAF_FAIL_Pop2+...)
make_typed_npass_tag() {
  local base="$1"      # e.g. MAF
  local thr="$2"       # threshold used to decide whether this tag is enabled
  local suffix="$3"    # e.g. SNP
  local pops_file="$4"

  local sum_expr
  sum_expr=$(
    cut -f2 "$pops_file" \
      | tr ',' '\n' \
      | sed '/^$/d' \
      | sort -u \
      | awk -v base="$base" -v suffix="$suffix" '
          BEGIN { first=1 }
          {
            if (!first) printf "+"
            printf "%s_PASS_%s_%s", base, suffix, $1
            first=0
          }'
  )

  printf 'NPASS_%s_%s:1=int(%s)\n' "$base" "$suffix" "$sum_expr"
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
  "$(make_typed_pass_tags "ExcHet" ">=" "${EH_POP_SNP:-}"        "SNP"       sample_groups.tsv)" \
  "$(make_typed_pass_tags "ExcHet" ">=" "${EH_POP_INDEL:-}"      "INDEL"     sample_groups.tsv)" \
  "$(make_typed_pass_tags "ExcHet" ">=" "${EH_POP_INVARIANT:-}"  "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_pass_tags "HWE"    ">=" "${HWE_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_pass_tags "HWE"    ">=" "${HWE_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_pass_tags "HWE"    ">=" "${HWE_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_pass_tags "MAF"    ">=" "${MAF_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_pass_tags "MAF"    ">=" "${MAF_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_pass_tags "MAF"    ">=" "${MAF_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_pass_tags "NS"     ">=" "${MIN_SAMPLES_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_pass_tags "NS"     ">=" "${MIN_SAMPLES_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_pass_tags "NS"     ">=" "${MIN_SAMPLES_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_pass_tags "CR"     ">=" "${MIN_CALLRATE_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_pass_tags "CR"     ">=" "${MIN_CALLRATE_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_pass_tags "CR"     ">=" "${MIN_CALLRATE_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)"
)

# Set up pass counting tags for each filter and concatenate together into single string for 1 pass annotation
# Only make npass tags if corrsponding pass_tags is created, i.e. filter is not disabled
# This avoids error in fill-tags
NPASS_TAGS=$(join_tags \
  "$(make_typed_npass_tag "ExcHet" "${EH_POP_SNP:-}"        "SNP"       sample_groups.tsv)" \
  "$(make_typed_npass_tag "ExcHet" "${EH_POP_INDEL:-}"      "INDEL"     sample_groups.tsv)" \
  "$(make_typed_npass_tag "ExcHet" "${EH_POP_INVARIANT:-}"  "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_npass_tag "HWE"    "${HWE_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_npass_tag "HWE"    "${HWE_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_npass_tag "HWE"    "${HWE_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_npass_tag "MAF"    "${MAF_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_npass_tag "MAF"    "${MAF_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_npass_tag "MAF"    "${MAF_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_npass_tag "NS"     "${MIN_SAMPLES_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_npass_tag "NS"     "${MIN_SAMPLES_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_npass_tag "NS"     "${MIN_SAMPLES_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)" \
  "$(make_typed_npass_tag "CR"     "${MIN_CALLRATE_POP_SNP:-}"       "SNP"       sample_groups.tsv)" \
  "$(make_typed_npass_tag "CR"     "${MIN_CALLRATE_POP_INDEL:-}"     "INDEL"     sample_groups.tsv)" \
  "$(make_typed_npass_tag "CR"     "${MIN_CALLRATE_POP_INVARIANT:-}" "INVARIANT" sample_groups.tsv)"
)

echo "POP_PASS_TAGS=$POP_PASS_TAGS"
echo "NPASS_TAGS=$NPASS_TAGS" 

# Subset to sites with just 2 alleles and samples above missing data filter
# Then add global annotations using fill-tags
# Then add per-pop threshold annotations using fill-tags
# Then add per-population pass flags, using the per-population threshold annotations using fill-tags
# Then count the number of populations that pass
# Note MAC is calculated from MAF (7 decimal precision), this could cause rounding for very large cohorts (i.e. 100k+)
# NOTE: bcftools +fill-tags breaks with any samples that have 2 letter names
bcftools view --threads ${1} -S ${4}.samples.txt -m2 -M2 -Ou "${3}" \
  | bcftools +setGT -Ou -- \
    -t q \
    -n . \
    -i "FORMAT/GQ < ${GENOTYPE_QUAL:-0} | FORMAT/DP < ${GENOTYPE_DP_MIN:-0} | FORMAT/DP > ${GENOTYPE_DP_MAX:-999999999}" \
  | bcftools +fill-tags -Ou - -- \
    -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING' \
  | bcftools +fill-tags -Ou - -- \
    -S sample_groups.tsv \
    -t 'NS,MAF,HWE,ExcHet,CR:1=1-F_MISSING' \
  | bcftools +fill-tags -Ou - -- \
    -t "$POP_PASS_TAGS" \
  | bcftools +fill-tags -Ou - -- \
    -t "$NPASS_TAGS" \
  | bcftools filter -Ou --SnpGap ${DIST_INDEL_GLOBAL_SNP:--999999999} --IndelGap ${DIST_INDEL_GLOBAL_INDEL:--999999999} \
  | bcftools filter -Ou -s MASK_FAIL -m+ -M vcf_masks.bed \
  | bcftools filter -Ou -s SNP_QUAL_FAIL   -m+ -e 'INFO/TYPE="SNP" && QUAL < '"${QUAL_GLOBAL_SNP:-0}" \
  | bcftools filter -Ou -s INDEL_QUAL_FAIL -m+ -e 'INFO/TYPE="INDEL" && QUAL < '"${QUAL_GLOBAL_INDEL:-0}" \
  | bcftools filter -Ou -s INV_QUAL_FAIL   -m+ -e 'INFO/TYPE="REF" && QUAL < '"${QUAL_GLOBAL_INVARIANT:-0}" \
  | bcftools filter -Ou -s SNP_DP_FAIL     -m+ -e 'INFO/TYPE="SNP" && (INFO/DP < '"${DP_MIN_GLOBAL_SNP:-0}"' || INFO/DP < '"${DP_LOWER_PERC_GLOBAL_SNP:-0}"' || INFO/DP > '"${DP_UPPER_PERC_GLOBAL_SNP:-999999999}"')' \
  | bcftools filter -Ou -s INDEL_DP_FAIL   -m+ -e 'INFO/TYPE="INDEL" && (INFO/DP < '"${DP_MIN_GLOBAL_INDEL:-0}"' || INFO/DP < '"${DP_LOWER_PERC_GLOBAL_INDEL:-0}"' || INFO/DP > '"${DP_UPPER_PERC_GLOBAL_INDEL:-999999999}"')' \
  | bcftools filter -Ou -s INV_DP_FAIL     -m+ -e 'INFO/TYPE="REF" && (INFO/DP < '"${DP_MIN_GLOBAL_INVARIANT:-0}"' || INFO/DP < '"${DP_LOWER_PERC_GLOBAL_INVARIANT:-0}"' || INFO/DP > '"${DP_UPPER_PERC_GLOBAL_INVARIANT:-999999999}"')' \
  | bcftools filter -Ou -s SNP_EH_FAIL     -m+ -e 'INFO/TYPE="SNP" && INFO/ExcHet < '"${EH_GLOBAL_SNP:--1}" \
  | bcftools filter -Ou -s INDEL_EH_FAIL   -m+ -e 'INFO/TYPE="INDEL" && INFO/ExcHet < '"${EH_GLOBAL_INDEL:--1}" \
  | bcftools filter -Ou -s SNP_HWE_FAIL    -m+ -e 'INFO/TYPE="SNP" && INFO/HWE < '"${HWE_GLOBAL_SNP:--1}" \
  | bcftools filter -Ou -s INDEL_HWE_FAIL  -m+ -e 'INFO/TYPE="INDEL" && INFO/HWE < '"${HWE_GLOBAL_INDEL:--1}" \
  | bcftools filter -Ou -s SNP_MAF_FAIL    -m+ -e 'INFO/TYPE="SNP" && INFO/MAF < '"${MAF_GLOBAL_SNP:-0}" \
  | bcftools filter -Ou -s INDEL_MAF_FAIL  -m+ -e 'INFO/TYPE="INDEL" && INFO/MAF < '"${MAF_GLOBAL_INDEL:-0}" \
  | bcftools filter -Ou -s SNP_NS_FAIL     -m+ -e 'INFO/TYPE="SNP" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_SNP:-0}" \
  | bcftools filter -Ou -s INDEL_NS_FAIL   -m+ -e 'INFO/TYPE="INDEL" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_INDEL:-0}" \
  | bcftools filter -Ou -s INV_NS_FAIL     -m+ -e 'INFO/TYPE="REF" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_INVARIANT:-0}" \
  | bcftools filter -Ou -s SNP_CR_FAIL     -m+ -e 'INFO/TYPE="SNP" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_SNP:-0}" \
  | bcftools filter -Ou -s INDEL_CR_FAIL   -m+ -e 'INFO/TYPE="INDEL" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_INDEL:-0}" \
  | bcftools filter -Ou -s INV_CR_FAIL     -m+ -e 'INFO/TYPE="REF" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_INVARIANT:-0}" \
  | bcftools filter -Ou -s SNP_POP_EH_FAIL   -m+ -e "${MIN_PASS:+INFO/TYPE=\"SNP\" && INFO/NPASS_ExcHet_SNP < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INDEL_POP_EH_FAIL -m+ -e "${MIN_PASS:+INFO/TYPE=\"INDEL\" && INFO/NPASS_ExcHet_INDEL < ${MIN_PASS}}" \
  | bcftools filter -Ou -s SNP_POP_HWE_FAIL   -m+ -e "${MIN_PASS:+INFO/TYPE=\"SNP\" && INFO/NPASS_HWE_SNP < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INDEL_POP_HWE_FAIL -m+ -e "${MIN_PASS:+INFO/TYPE=\"INDEL\" && INFO/NPASS_HWE_INDEL < ${MIN_PASS}}" \
  | bcftools filter -Ou -s SNP_POP_MAF_FAIL   -m+ -e "${MIN_PASS:+INFO/TYPE=\"SNP\" && INFO/NPASS_MAF_SNP < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INDEL_POP_MAF_FAIL -m+ -e "${MIN_PASS:+INFO/TYPE=\"INDEL\" && INFO/NPASS_MAF_INDEL < ${MIN_PASS}}" \
  | bcftools filter -Ou -s SNP_POP_NS_FAIL    -m+ -e "${MIN_PASS:+INFO/TYPE=\"SNP\" && INFO/NPASS_NS_SNP < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INDEL_POP_NS_FAIL  -m+ -e "${MIN_PASS:+INFO/TYPE=\"INDEL\" && INFO/NPASS_NS_INDEL < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INV_POP_NS_FAIL    -m+ -e "${MIN_PASS:+INFO/TYPE=\"REF\" && INFO/NPASS_NS_INVARIANT < ${MIN_PASS}}" \
  | bcftools filter -Ou -s SNP_POP_CR_FAIL    -m+ -e "${MIN_PASS:+INFO/TYPE=\"SNP\" && INFO/NPASS_CR_SNP < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INDEL_POP_CR_FAIL  -m+ -e "${MIN_PASS:+INFO/TYPE=\"INDEL\" && INFO/NPASS_CR_INDEL < ${MIN_PASS}}" \
  | bcftools filter -Ou -s INV_POP_CR_FAIL    -m+ -e "${MIN_PASS:+INFO/TYPE=\"REF\" && INFO/NPASS_CR_INVARIANT < ${MIN_PASS}}" \
  | bcftools view --threads ${1} -Ob -o tmp.bcf

# Drop failing sites to create filtered vcf file (main output)
bcftools view --threads "${1}" -f PASS -Ou tmp.bcf \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' -Oz9 -o ${4}.filt.vcf.gz 
bcftools index --threads ${1} -t ${4}.filt.vcf.gz

# Create filtered sitelist file by dropping non-passing variants and genotypes
bcftools view --threads "${1}" -G -f PASS -Ou tmp.bcf \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' -Oz9 -o ${4}.sitelist.vcf.gz 
bcftools index --threads ${1} -t ${4}.sitelist.vcf.gz

# Create metrics file for QC histograms
mapfile -t INFO_TAGS < <(
  bcftools view -h tmp.bcf \
    | awk -F'[=,]' '/^##INFO=<ID=/{print $3}' \
    | awk '
        $0 == "TYPE" ||
        $0 == "DP" ||
        $0 == "ExcHet" ||
        $0 == "HWE" ||
        $0 == "MAF" ||
        $0 == "NS" ||
        $0 == "CR" ||
        $0 ~ /^(NS|MAF|HWE|ExcHet|CR)_[^_]+$/'
)

HEADER="CHROM\tPOS\tFILTER\tQUAL"
FMT='%CHROM\t%POS\t%FILTER\t%QUAL'

for tag in "${INFO_TAGS[@]}"; do
  HEADER="${HEADER}\t${tag}"
  FMT="${FMT}\t%INFO/${tag}"
done
FMT="${FMT}\n"

{
  printf '%b\n' "$HEADER"
  bcftools query -f "$FMT" tmp.bcf
} | bgzip -c > "${4}.metrics.tsv.gz"

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n ${4}.filt.vcf.gz | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${4}.counts"