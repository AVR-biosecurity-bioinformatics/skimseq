#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = interval_hash

base="$(basename ${3})"
prefix="${base%.vcf.gz}"

# Create a small summary of the number of sites passing and failing each filter
bcftools query -f '%FILTER\n' ${3} \
  | sort \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"} {print $2, $1}' \
  > "${prefix}_filter_summary.tsv"

# ------- make filter summary histograms ------

# We bin the values and create the histogram in awk to avoid parsing massive files to R

# ---- helper functions ----
# Compute NBINS-bin width & origin (SITE metric format like "%INFO/QD\n")
fixed_bins_site() {
  local INPUT="$1" FMT="$2" NBINS="${3:-50}"
  (( NBINS < 1 )) && NBINS=1

  # Emit *something* even if bcftools fails
  local tmp out rc=0
  tmp=$(mktemp) || { echo "1 0"; return; }

  # Collect values; ignore missing; compute in awk
  if ! bcftools query -f "$FMT" "$INPUT" 2>/dev/null \
      | awk -v NB="$NBINS" '
          $1!="." && $1!="" { v=$1+0; n++; if(n==1){min=v;max=v}else{if(v<min)min=v;if(v>max)max=v} }
          END{
            if(n<2 || max<=min){ print 1, (n?min:0); exit }
            bw=(max-min)/NB; if(bw<=0)bw=1
            printf("%.12g %.12g\n", bw, min)
          }' > "$tmp"; then
    rc=1
  fi

  out=$(cat "$tmp"); rm -f "$tmp"
  if [[ $rc -ne 0 || -z $out ]]; then
    echo "1 0"    # safe fallback: width=1, origin=0
  else
    echo "$out"
  fi
}

# Compute NBINS-bin width & origin (GT metric tag like "%DP" / "%GQ")
fixed_bins_gt() {
  local INPUT="$1" METRIC="$2" NBINS="${3:-50}"
  (( NBINS < 1 )) && NBINS=1

  local tmp out rc=0
  tmp=$(mktemp) || { echo "1 0"; return; }

  if ! bcftools query -f "[${METRIC}\t]\n" "$INPUT" 2>/dev/null \
      | tr '\t' '\n' \
      | awk -v NB="$NBINS" '
          $1!="." && $1!="" { v=$1+0; n++; if(n==1){min=v;max=v}else{if(v<min)min=v;if(v>max)max=v} }
          END{
            if(n<2 || max<=min){ print 1, (n?min:0); exit }
            bw=(max-min)/NB; if(bw<=0)bw=1
            printf("%.12g %.12g\n", bw, min)
          }' > "$tmp"; then
    rc=1
  fi

  out=$(cat "$tmp"); rm -f "$tmp"
  if [[ $rc -ne 0 || -z $out ]]; then
    echo "1 0"
  else
    echo "$out"
  fi
}

# Bin numeric stream using BIN width anchored at ORIGIN
bin_stream() {
  local BIN="$1" ORG="${2:-0}"
  awk -v BIN="$BIN" -v ORG="$ORG" '
    $1=="." || $1=="" { next }
    { v=$1+0; b = int((v-ORG)/BIN)*BIN + ORG; print b }
  '
}

# Collapse bins to "RULE  FILTER  VTYPE  BIN  COUNT"
emit_counts() {
  awk -v rule="$1" -v status="$2" -v vt="$3" -v OFS='\t' '
    { h[$1]++ }
    END { n=0; for (b in h) bins[n++]=b; asort(bins);
          for(i=1;i<=n;i++){ b=bins[i]; print rule, status, vt, b+0, h[b] } }'
}

# create_pf_histogram MODE INPUT FAILSEL METRIC RULELABEL VTYPE [NBINS]
create_pf_histogram() {
  local MODE="$1" INPUT="$2" FAILSEL="$3" METRIC="$4" RULE="$5" VTYPE="$6" NBINS="${7:-50}"

  if [[ "$MODE" == "SITE" ]]; then
    local BW ORG
    # Protect against empty output from helper under `set -u`
    read -r BW ORG < <(fixed_bins_site "$INPUT" "$METRIC" "$NBINS")
    BW=${BW:-1}; ORG=${ORG:-0}

    bcftools query -i "FILTER~\"${FAILSEL}\"" -f "$METRIC" "$INPUT" \
      | bin_stream "$BW" "$ORG" \
      | emit_counts "$RULE" "FAIL" "$VTYPE"

    bcftools query -i "FILTER!~\"${FAILSEL}\"" -f "$METRIC" "$INPUT" \
      | bin_stream "$BW" "$ORG" \
      | emit_counts "$RULE" "PASS" "$VTYPE"

  elif [[ "$MODE" == "GT" ]]; then
    local BW ORG
    read -r BW ORG < <(fixed_bins_gt "$INPUT" "$METRIC" "$NBINS")
    BW=${BW:-1}; ORG=${ORG:-0}

    bcftools query -f "[${METRIC},%FT\t]\n" "$INPUT" \
      | tr '\t' '\n' \
      | awk -F',' -v p="$FAILSEL" '$1!="." && $1!="" && $2 ~ p { print $1 }' \
      | bin_stream "$BW" "$ORG" \
      | emit_counts "$RULE" "FAIL" "$VTYPE"

    bcftools query -f "[${METRIC},%FT\t]\n" "$INPUT" \
      | tr '\t' '\n' \
      | awk -F',' -v p="$FAILSEL" '$1!="." && $1!="" && $2 !~ p { print $1 }' \
      | bin_stream "$BW" "$ORG" \
      | emit_counts "$RULE" "PASS" "$VTYPE"
  else
    echo "create_pf_histogram: unknown MODE '$MODE'" >&2
    return 2
  fi
}

# ---- build the table ----
out="${prefix}_filter_hist.tsv"
printf "RULE\tFILTER\tVARIANT_TYPE\tBIN\tCOUNT\n" > "$out"

VTYPE="${4}"  # snp|indel|invariant
NBINS=100 # Maximum number of data bins

# Site-level histograms (use tmp.tagged.bcf)
INPUT_SITE=${3}

# Global filters
create_pf_histogram SITE "$INPUT_SITE" "QUAL_FAIL"  "%QUAL\n"                QUAL        "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "EH_FAIL"    "%INFO/ExcHet\n"      ExcHet   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "HWE_FAIL"    "%INFO/HWE\n"      HWE   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "DP_FAIL"    "%INFO/DP\n"             DP          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "DIST_INDEL_FAIL"    "%INFO/DIST_INDEL\n"      DIST_INDEL          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "MAF_FAIL"   "%INFO/MAF\n"            MAF         "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "NS_FAIL"    "%INFO/NS\n"             NS          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "CR_FAIL"    "%INFO/CR\n"             CR          "$VTYPE" "$NBINS" >> "$out"

# pop filters
create_pf_histogram SITE "$INPUT_SITE" "POP_EH_FAIL"    "%INFO/ExcHet\n"      ExcHet   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "POP_HWE_FAIL"    "%INFO/HWE\n"      HWE_POP   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "POP_MAF_FAIL"   "%INFO/MAF\n"            MAF_POP         "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "POP_NS_FAIL"    "%INFO/NS\n"             NS_POP          "$VTYPE" "$NBINS" >> "$out"

# Zip output summary table
pigz -p ${1} $out

