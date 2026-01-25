#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed
# $6 = interval_hash
# $7 = DP summary

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

# Calculate percentile DP filters from DP histogram
read DPlower DPupper < <(
  awk -v pl="$PCT_LOW" -v ph="$PCT_HIGH" '
    { dp[NR]=$1; cnt[NR]=$2+0; N+=cnt[NR] }
    END{
      if(N==0) exit 1
      low  = pl/100*N; li=int(low);  if(li<low)  li++; if(li<1) li=1; if(li>N) li=N
      high = ph/100*N; ui=int(high); if(ui<high) ui++; if(ui<1) ui=1; if(ui>N) ui=N
      cum=0
      for(i=1;i<=NR;i++){
        cum += cnt[i]
        if(!lo && cum>=li) lo=dp[i]
        if(!hi && cum>=ui){ hi=dp[i]; break }
      }
      printf "%d %d\n", lo, hi
    }
  ' "$7"
)

# Subset to target variant class and run drop all genotypes, then run site-level soft filtering 
# (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)
bcftools view --threads ${1} -G ${TYPE_ARGS} -Ou "${3}" \
  | bcftools filter -Ou -s QUAL_FAIL   -m+ -e "QUAL     <= ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s EH_FAIL     -m+ -e "INFO/ExcHet <= ${EH:-1e9}" \
  | bcftools filter -Ou -s HWE_FAIL     -m+ -e "INFO/HWE <= ${HWE:-1e9}" \
  | bcftools filter -Ou -s DP_FAIL     -m+ -e "INFO/DP <= ${DPmin:-0} || INFO/DP <= ${DPlower:-0} || INFO/DP >= ${DPupper:-999999999}" \
  | bcftools filter -Ou -s INDELDIST_FAIL   -m+ -e "INFO/DIST_INDEL <= ${DIST_INDEL:-0}" \
  | bcftools filter -Ou -s MAF_FAIL    -m+ -e "INFO/MAF <= ${MAF:-0}" \
  | bcftools filter -Ou -s MAC_FAIL    -m+ -e "INFO/MAC <= ${MAC:-0}" \
  | bcftools filter -Ou -s NS_FAIL     -m+ -e "INFO/NS <= ${NS:-0}" \
  | bcftools filter -Ou -s CR_FAIL     -m+ -e "INFO/CR <= ${CR:-0}" \
  | bcftools filter -Ou -s MASK_FAIL   -m+ -M vcf_masks.bed \
  | bcftools view --threads ${1} -Ob -o tmp.tagged.bcf

# Keep only variants that PASS & index output
# TODO: Drop FT and other extra fields from vcf
bcftools view --threads ${1} -f PASS -Oz -o ${6}_${4}_filtered.vcf.gz tmp.tagged.bcf
bcftools index --threads ${1} -t ${6}_${4}_filtered.vcf.gz

# Output number of variant records remaining (non-header lines)
nvars=$(bcftools index -n "${6}_${4}_filtered.vcf.gz" | tr -d '[:space:]')
printf "%s\n" "$nvars" > "${6}_${4}.counts"

# Create a small summary of the number of sites passing and failing each filter
bcftools query -f '%FILTER\n' tmp.tagged.bcf \
  | sort \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"} {print $2, $1}' \
  > "${6}_${4}_filter_summary.tsv"

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
out="${6}_${4}_filter_hist.tsv"
printf "RULE\tFILTER\tVARIANT_TYPE\tBIN\tCOUNT\n" > "$out"

VTYPE="${4}"  # snp|indel|invariant
NBINS=100 # Maximum number of data bins

# Site-level histograms (use tmp.tagged.bcf)
INPUT_SITE=tmp.tagged.bcf
create_pf_histogram SITE "$INPUT_SITE" "QUAL_FAIL"  "%QUAL\n"                QUAL        "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "EH_FAIL"    "%INFO/ExcHet\n"      ExcHet   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "HWE_FAIL"    "%INFO/HWE\n"      HWE   "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "DP_FAIL"    "%INFO/DP\n"             DP          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "DIST_INDEL_FAIL"    "%INFO/DIST_INDEL\n"      DIST_INDEL          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "MAF_FAIL"   "%INFO/MAF\n"            MAF         "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "MAC_FAIL"   "%INFO/MAC\n"            MAC         "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "NS_FAIL"    "%INFO/NS\n"             NS          "$VTYPE" "$NBINS" >> "$out"
create_pf_histogram SITE "$INPUT_SITE" "CR_FAIL"    "%INFO/CR\n"             CR          "$VTYPE" "$NBINS" >> "$out"

# Zip output summary table
pigz -p ${1} $out

# Remove temporary vcf files
rm -f tmp* *.bcf