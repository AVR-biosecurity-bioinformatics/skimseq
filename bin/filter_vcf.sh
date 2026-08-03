#!/bin/bash
set -uoe pipefail

# Function to check pipeline
check_pipeline() {
    local -a st=("$@")

    echo "PIPESTATUS: ${st[*]}" >&2

    for i in "${!st[@]}"; do
        if (( st[i] != 0 )); then
            echo "Pipeline stage $((i + 1)) failed with exit code ${st[i]}" >&2
            return "${st[i]}"
        fi
    done

    return 0
}

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

mapfile -t POPS < <(
  cut -f2 sample_groups.tsv \
    | tr ',' '\n' \
    | sed '/^$/d' \
    | sort -u
)

N_POPS=${#POPS[@]}

if [[ "$N_POPS" -eq 0 ]]; then
  echo "ERROR: no populations retained in sample_groups.tsv after filtering" >&2
  exit 1
fi


# Population-level filtering mode:
#   ALL = fail site only if all populations fail threshold
#   ANY = fail site if any population fails threshold
POPULATION_FAIL_MODE="${POPULATION_FAIL_MODE:-ALL}"

case "$POPULATION_FAIL_MODE" in
  ALL|ANY) ;;
  *)
    echo "ERROR: POPULATION_FAIL_MODE must be ALL or ANY, got '$POPULATION_FAIL_MODE'" >&2
    exit 1
    ;;
esac

# Function to join non-empty strings with a separator
join_by() {
  local sep="$1"
  shift
  local out=""
  local x
  for x in "$@"; do
    [[ -z "$x" ]] && continue
    [[ -n "$out" ]] && out+="$sep"
    out+="$x"
  done
  printf '%s\n' "$out"
}

# Function to Return a per-type population clause, e.g.
#   INFO/TYPE="SNP" && INFO/MAF_Pop1<0.05 && INFO/MAF_Pop2<0.05
# or
#   INFO/TYPE="SNP" && (INFO/MAF_Pop1<0.05 || INFO/MAF_Pop2<0.05)
#
# If threshold is disabled, return empty string.
make_type_pop_clause() {
  local base="$1"   # MAF / HWE / ExcHet / NS / CR
  local op="$2"     # usually <
  local thr="$3"
  local type="$4"   # SNP / INDEL / REF

  if [[ -z "${thr:-}" || "$thr" == "NA" || "$thr" == "na" ]]; then
    printf ''
    return 0
  fi

  local parts=()
  local pop
  for pop in "${POPS[@]}"; do
    parts+=( "INFO/${base}_${pop}${op}${thr}" )
  done

  local expr
  if [[ "$POPULATION_FAIL_MODE" == "ALL" ]]; then
    expr=$(join_by ' && ' "${parts[@]}")
    printf 'INFO/TYPE="%s" && %s' "$type" "$expr"
  else
    expr=$(join_by ' || ' "${parts[@]}")
    printf 'INFO/TYPE="%s" && (%s)' "$type" "$expr"
  fi
}

# Build the final metric-level population expression by combining type clauses.
# Example:
#   make_metric_pop_expr "MAF" "<" "SNP:0.05" "INDEL:0.05"
#
# returns something like:
#   (INFO/TYPE="SNP" && (...)) || (INFO/TYPE="INDEL" && (...))
#
# If all thresholds are disabled, returns an impossible false expression.
make_metric_pop_expr() {
  local base="$1"
  local op="$2"
  shift 2

  local clauses=()
  local spec type thr clause

  for spec in "$@"; do
    type="${spec%%:*}"
    thr="${spec#*:}"
    clause=$(make_type_pop_clause "$base" "$op" "$thr" "$type")
    [[ -n "$clause" ]] && clauses+=( "( ${clause} )" )
  done

  if [[ "${#clauses[@]}" -eq 0 ]]; then
    # valid no-op false expression
    printf 'INFO/TYPE="SNP" && INFO/TYPE="INDEL"'
    return 0
  fi

  printf '%s\n' "$(join_by ' || ' "${clauses[@]}")"
}

# Build global filtering expressions
QUAL_EXPR='(INFO/TYPE="SNP" && QUAL < '"${QUAL_GLOBAL_SNP:-0}"') || (INFO/TYPE="INDEL" && QUAL < '"${QUAL_GLOBAL_INDEL:-0}"') || (INFO/TYPE="REF" && QUAL < '"${QUAL_GLOBAL_INVARIANT:-0}"')'

DP_MIN_EXPR='( INFO/TYPE="SNP" && INFO/DP < '"${DP_MIN_GLOBAL_SNP:-0}"') || ( INFO/TYPE="INDEL" && INFO/DP < '"${DP_MIN_GLOBAL_INDEL:-0}"') || (INFO/TYPE="REF"   && INFO/DP < '"${DP_MIN_GLOBAL_INVARIANT:-0}"')'

DP_LOWER_PERC_EXPR='(INFO/TYPE="SNP" && INFO/DP < '"${DP_LOWER_PERC_GLOBAL_SNP:-0}"') || (INFO/TYPE="INDEL" && INFO/DP < '"${DP_LOWER_PERC_GLOBAL_INDEL:-0}"') || ( INFO/TYPE="REF" && INFO/DP < '"${DP_LOWER_PERC_GLOBAL_INVARIANT:-0}"')'

DP_UPPER_PERC_EXPR='(INFO/TYPE="SNP" && INFO/DP > '"${DP_UPPER_PERC_GLOBAL_SNP:-999999999}"') || (INFO/TYPE="INDEL" && INFO/DP > '"${DP_UPPER_PERC_GLOBAL_INDEL:-999999999}"') || (INFO/TYPE="REF" && INFO/DP > '"${DP_UPPER_PERC_GLOBAL_INVARIANT:-999999999}"')'

EH_EXPR='(INFO/TYPE="SNP" && INFO/ExcHet < '"${EH_GLOBAL_SNP:--1}"') || (INFO/TYPE="INDEL" && INFO/ExcHet < '"${EH_GLOBAL_INDEL:--1}"')' 

HWE_EXPR='(INFO/TYPE="SNP" && INFO/HWE < '"${HWE_GLOBAL_SNP:--1}"') || (INFO/TYPE="INDEL" && INFO/HWE < '"${HWE_GLOBAL_INDEL:--1}"')'

MAF_EXPR='(INFO/TYPE="SNP" && INFO/MAF < '"${MAF_GLOBAL_SNP:-0}"') || (INFO/TYPE="INDEL" && INFO/MAF < '"${MAF_GLOBAL_INDEL:-0}"')'

NS_EXPR='(INFO/TYPE="SNP" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_SNP:-0}"') || (INFO/TYPE="INDEL" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_INDEL:-0}"') || (INFO/TYPE="REF" && INFO/NS < '"${MIN_SAMPLES_GLOBAL_INVARIANT:-0}"')'

CR_EXPR='(INFO/TYPE="SNP" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_SNP:-0}"') || (INFO/TYPE="INDEL" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_INDEL:-0}"') || (INFO/TYPE="REF" && INFO/CR < '"${MIN_CALLRATE_GLOBAL_INVARIANT:-0}"')'

# Build per-pop filtering expressions
POP_EH_EXPR=$(make_metric_pop_expr "ExcHet" "<" \
  "SNP:${EH_POP_SNP:-}" \
  "INDEL:${EH_POP_INDEL:-}")

POP_HWE_EXPR=$(make_metric_pop_expr "HWE" "<" \
  "SNP:${HWE_POP_SNP:-}" \
  "INDEL:${HWE_POP_INDEL:-}")

POP_MAF_EXPR=$(make_metric_pop_expr "MAF" "<" \
  "SNP:${MAF_POP_SNP:-}" \
  "INDEL:${MAF_POP_INDEL:-}")

POP_NS_EXPR=$(make_metric_pop_expr "NS" "<" \
  "SNP:${MIN_SAMPLES_POP_SNP:-}" \
  "INDEL:${MIN_SAMPLES_POP_INDEL:-}" \
  "REF:${MIN_SAMPLES_POP_INVARIANT:-}")

POP_CR_EXPR=$(make_metric_pop_expr "CR" "<" \
  "SNP:${MIN_CALLRATE_POP_SNP:-}" \
  "INDEL:${MIN_CALLRATE_POP_INDEL:-}" \
  "REF:${MIN_CALLRATE_POP_INVARIANT:-}")

# Subset to sites with just 2 alleles and samples above missing data filter
# Then add global annotations using fill-tags
# Then add per-pop threshold annotations using fill-tags
# Note MAC is calculated from MAF (7 decimal precision), this could cause rounding for very large cohorts (i.e. 100k+)
# NOTE: bcftools +fill-tags breaks with any samples that have 2 letter names
set +e
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
  | bcftools filter -Ou --SnpGap ${DIST_INDEL_GLOBAL_SNP:--999999999} --IndelGap ${DIST_INDEL_GLOBAL_INDEL:--999999999} \
  | bcftools filter -Ou -s MASK_FAIL -m+ -M vcf_masks.bed \
  | bcftools filter -Ou -s QUAL_FAIL -m+ -e "$QUAL_EXPR" \
  | bcftools filter -Ou -s DP_MIN_FAIL -m+ -e "$DP_MIN_EXPR" \
  | bcftools filter -Ou -s DP_LOWER_PERC_FAIL -m+ -e "$DP_LOWER_PERC_EXPR" \
  | bcftools filter -Ou -s DP_UPPER_PERC_FAIL -m+ -e "$DP_UPPER_PERC_EXPR" \
  | bcftools filter -Ou -s EH_FAIL -m+ -e "$EH_EXPR" \
  | bcftools filter -Ou -s HWE_FAIL -m+ -e "$HWE_EXPR" \
  | bcftools filter -Ou -s MAF_FAIL -m+ -e "$MAF_EXPR" \
  | bcftools filter -Ou -s NS_FAIL -m+ -e "$NS_EXPR" \
  | bcftools filter -Ou -s CR_FAIL -m+ -e "$CR_EXPR" \
  | bcftools filter -Ou -s POP_EH_FAIL  -m+ -e "$POP_EH_EXPR" \
  | bcftools filter -Ou -s POP_HWE_FAIL -m+ -e "$POP_HWE_EXPR" \
  | bcftools filter -Ou -s POP_MAF_FAIL -m+ -e "$POP_MAF_EXPR" \
  | bcftools filter -Ou -s POP_NS_FAIL  -m+ -e "$POP_NS_EXPR" \
  | bcftools filter -Ou -s POP_CR_FAIL  -m+ -e "$POP_CR_EXPR" \
  | bcftools view --threads ${1} -Ob -o tmp.bcf

# Catch error codes from piped tools so nextflow can retry
st=("\${PIPESTATUS[@]}")
set -e
check_pipeline "\${st[@]}" || exit \$?

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