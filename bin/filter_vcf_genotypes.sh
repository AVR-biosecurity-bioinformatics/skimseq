#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf

# Add FORMAT/FT tags using awk and annotate - BCFtools doesnt natively support soft filtering of genotypes
bcftools query -f '%CHROM\t%POS[\t%GQ\t%DP]\n' ${3} \
| awk -v OFS='\t' -v gq="${GQ:-0}" -v dmin="${gtDPmin:-0}" -v dmax="${gtDPmax:-999999999}" '
  {
    printf "%s%s%s", $1, OFS, $2
    # fields: 3=GQ_s1, 4=DP_s1, 5=GQ_s2, 6=DP_s2, ...
    for (i=3; i<=NF; i+=2) {
      g=$i; d=$(i+1)
      if (g=="." || d==".") {
        ft="."
      } else {
        ft=""
        if (g+0 < gq)  ft = (ft=="" ? "GQ_FAIL"   : ft ";GQ_FAIL")
        if (d+0 < dmin)ft = (ft=="" ? "GTDP_FAIL": ft ";GTDP_FAIL")
        if (d+0 > dmax)ft = (ft=="" ? "GTDP_FAIL": ft ";GTDP_FAIL")
        if (ft=="") ft="PASS"
      }
      printf OFS "%s", ft
    }
    print ""
  }' \
| bgzip -c > FT.tsv.gz

tabix -s1 -b2 -e2 FT.tsv.gz

# add header + inject FORMAT/FT
cat > FT.hdr <<'EOF'
##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filters (per-sample) from GQ/DP thresholds">
EOF

# Annotate filter column in vcf
bcftools annotate \
  -h FT.hdr \
  -a FT.tsv.gz \
  -c CHROM,POS,FORMAT/FT \
  -Ob -o gt_masked.bcf ${3}

# Set failing GTs to missing 
bcftools +setGT -Ou gt_masked.bcf -- -t q -n . -e 'FMT/FT="PASS" && FMT/FT="."' \
  | bcftools view -U -S samples_to_keep.txt -Ob gt_filtered.bcf

# Calculate per_sample missing data
bcftools stats -s - gt_filtered.bcf \
| awk -v out="missing_summary.tsv" 'BEGIN{OFS="\t"}
    $1=="SN" && $3=="number" && $5=="records:" {
        total=$6
        next
    }
    $1=="PSC" {
        if(!printed_header){
            print "#TOTAL_RECORDS", total > out
            print "SAMPLE","NMISS" > out
            printed_header=1
        }
        print $3, $14 > out
    }
'
# Find samples above the missing fraction filter
awk -v thr="$MISSING_FRAC" 'NR==1 {next} $4!="NA" && ($4+0) < thr {print $1}' \
missing_summary.tsv > samples_to_keep.txt
  
# Re-calculate minor alelle count (MAC) info tag
# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp.bcf \
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

# Annotate the vcf with INFO/MAC and update site tags
bcftools annotate -h MAC.hdr -a MAC.tsv.gz -c CHROM,POS,INFO/MAC -Ou tmp.bcf  \
  | bcftools +fill-tags -- -t MAF,ExcHet,HWE,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \
  | bcftools view -U -Oz9 -o final.vcf

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

# Genotype-level histograms (use gt_masked.bcf which has FORMAT/FT)
INPUT_GT=gt_masked.bcf
create_pf_histogram GT "$INPUT_GT" '(^|;)GQ_FAIL(;|$)'        "%GQ" GT_GQ "ALL" >> "$out"
create_pf_histogram GT "$INPUT_GT" '(^|;)GTDP_FAIL(;|$)'      "%DP" GT_DP "ALL" >> "$out"


# Zip output summary table
pigz -p ${1} $out

# Remove temporary vcf files
rm -f tmp* MAC.tsv.gz* *.hdr *.bcf