#!/bin/bash
set -euo pipefail


## args
# $1 = cpus
# $2 = memory    (unused placeholder)
# $3 = sample
# $4 = interval_hash
# $5 = stderr.log (HaplotypeCaller)
# $6 = assembly.tsv (IGV-style activity track)
# $7 = BAM/CRAM (indexed)
# $8 = VCF/BCF   (bgz+index)

cpus="${1:-1}"
sample="$3"
ihash="$4"
logfile="$5"
assembly_tsv="$6"
bam="$7"
vcf="$8"

# Optional: set dedup handling to mirror NotDuplicateReadFilter (default: true)
dedup="${DEDUP:-true}"   # export DEDUP=false to allow duplicates
[[ "$dedup" == "false" ]] && SAMTOOLS_FFLAG=260 || SAMTOOLS_FFLAG=1284

# Tidy temp workspace
tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

###############################################################################
# 1) Progress → progress_summary.tsv  (one tiny file, easy to debug)
###############################################################################
printf "sample\tinterval_hash\tchrom\tpos\telapsed_min\tregions\tregions_per_min\n" > "$tmpdir/progress.tsv"
grep -F 'ProgressMeter -' "$logfile" \
| awk -v s="$sample" -v h="$ihash" '
  match($0,/ProgressMeter -[[:space:]]*([[:alnum:]_.]+):([0-9]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9.]+)/,m){
    printf("%s\t%s\t%s\t%s\t%.3f\t%d\t%.1f\n", s,h,m[1],m[2],m[3],m[4],m[5])
  }' >> "$tmpdir/progress.tsv"

###############################################################################
# 2) Assembly TSV → hc_ar.bed   (chrom start end active size)
###############################################################################
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 || $1 ~ /^#/ { next }
  $4=="end-marker"   { next }
  $4 ~ /^size=/ {
    sz=$4; sub(/^size=/,"",sz)
    act = ($5+0>0)?1:0
    print $1,$2,$3,act,sz
  }' "$assembly_tsv" \
| sort -k1,1 -k2,2n > "$tmpdir/hc_ar.bed"

###############################################################################
# 3) Tick windows (first window clamped to first region start on contig)
###############################################################################
awk -F'\t' '!seen[$1]++{print $1"\t"$2}' "$tmpdir/hc_ar.bed" > "$tmpdir/chrom_first_start.tsv"

awk -F'\t' -v OFS='\t' '
NR==FNR { firstS[$1]=$2; next }      # load first start per chrom
NR==1 { next }                       # skip header
{
  s=$1; h=$2; c=$3; p=$4+0; em=$5+0; r=$6+0;
  k=s "|" h "|" c
  if (!(k in seen)) {
    ws = (c in firstS ? firstS[c] : 0)
    ms = p-1; if (ms<0) ms=0
    if (ws>ms) ws=ms
    we = (p>ws? p : ws+1)
    dem=em; dr=r; idx[k]=1
    print c,ws,we,s,h,dem,dr,(s "|" h "|" c "|" idx[k])
  } else {
    ws = prevP[k]-1; if (ws<0) ws=0
    we = (p>ws? p : ws+1)
    dem = em-prevE[k]; if (dem<0) dem=0
    dr  = r -prevR[k]; if (dr <0) dr =0
    idx[k]++
    print c,ws,we,s,h,dem,dr,(s "|" h "|" c "|" idx[k])
  }
  seen[k]=1; prevP[k]=p; prevE[k]=em; prevR[k]=r
}' "$tmpdir/chrom_first_start.tsv" "$tmpdir/progress.tsv" > "$tmpdir/tick_windows.bed"
# cols: chrom start end sample interval_hash d_elapsed_min d_regions window_id

###############################################################################
# 4) Tag each assembly region with its finishing window
###############################################################################
bedtools intersect \
  -a <(awk 'BEGIN{OFS="\t"} {print $1,$3-1,$3,"RID"NR,$2,$3,$4,$5}' "$tmpdir/hc_ar.bed") \
  -b "$tmpdir/tick_windows.bed" -wao \
| awk -F'\t' -v OFS='\t' '
BEGIN{
  print "sample","interval_hash","window_id","window_start","window_end","window_elapsed_min","window_regions",
        "region_id","chrom","start","end","active","size_bp"
}
# A: 1 chromA 2 end-1 3 end 4 region_id 5 start 6 end 7 active 8 size
# B: 9 chromW 10 ws 11 we 12 sample 13 hash 14 d_elapsed 15 d_regions 16 window_id 17 overlap
$17+0>0 {
  print $12,$13,$16,$10,$11,$14,$15, $4,$1,$5,$6,$7,$8
}' > "$tmpdir/regions_with_windows.tsv"

# A compact BED for counting; keep "active"
awk -F'\t' -v OFS='\t' '
NR==1{next}
{ print $9,$10,$11,$8,$3,$1,$2,$4,$5,$6,$12 }' "$tmpdir/regions_with_windows.tsv" > "$tmpdir/regions_with_windows.bed"
# cols: 1 chrom 2 start 3 end 4 region_id 5 window_id 6 sample 7 interval_hash 8 window_start 9 window_end 10 window_elapsed_min 11 active

###############################################################################
# 5) VCF covariates (SNP/INDEL counts per region) → window summary
###############################################################################
bedtools intersect -a "$tmpdir/regions_with_windows.bed" \
  -b <(bcftools view -H -v snps   "$vcf" | awk -v OFS="\t" '{print $1,$2-1,$2}') -c \
| bedtools intersect -a stdin \
  -b <(bcftools view -H -v indels "$vcf" | awk -v OFS="\t" '{print $1,$2-1,$2}') -c \
| awk -F'\t' -v OFS='\t' '
BEGIN{
  print "sample","interval_hash","window_id","window_start","window_end","window_elapsed_min",
        "n_regions","n_active_regions","n_inactive_regions",
        "sum_variants","sum_indels","mean_variants_per_region","mean_indels_per_region"
}
{
  # fields after both -c intersects:
  # 1 chrom 2 start 3 end 4 region_id 5 window_id 6 sample 7 interval_hash 8 win_start 9 win_end 10 win_min 11 active 12 n_snps 13 n_indels
  wid=$5; s=$6; h=$7; wst=$8; wend=$9; wmin=$10
  act=($11==1); v=$12+0; ins=$13+0
  n[wid]++; a[wid]+=act; varsum[wid]+=v; indsum[wid]+=ins
  S[wid]=s; H[wid]=h; WST[wid]=wst; WEN[wid]=wend; WMIN[wid]=wmin
}
END{
  for (k in n){
    ninact = n[k]-a[k]
    mv=(n[k]? varsum[k]/n[k]:0)
    mi=(n[k]? indsum[k]/n[k]:0)
    print S[k],H[k],k,WST[k],WEN[k],WMIN[k],
          n[k],a[k],ninact,
          varsum[k],indsum[k],mv,mi
  }
}' > "$tmpdir/window_summary.tsv"

###############################################################################
# 6) BAM covariates per window (genomic_bases, aligned_bases, I/D/S)
###############################################################################
# windows list: chrom start end window_id sample hash elapsed_min regions
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8,$4,$5,$6,$7}' "$tmpdir/tick_windows.bed" > "$tmpdir/windows.list"

# genomic_bases from assembly regions (no per-window loop needed here)
# bedtools coverage: A=windows, B=hc_ar.bed → bases of A covered by B (clipped)
bedtools coverage -a "$tmpdir/tick_windows.bed" -b "$tmpdir/hc_ar.bed" \
| awk -F'\t' -v OFS='\t' '{print $8, $13}' > "$tmpdir/win_genomic_bases.tsv"
# columns: window_id  bp_covered_by_regions

# aligned_bases: sum(depth) over (regions ∩ windows). We do a small loop for I/D/S anyway;
# keep aligned_bases in that same loop for simplicity & correctness.
echo -e "window_id\tgenomic_bases\taligned_bases\tI_sum\tD_sum\tS_sum\tIDS_sum" > "$tmpdir/window_cov.tsv"

while IFS=$'\t' read -r chrom ws we wid sample ihash wmins wregs; do
  # Build union of (regions in this window) clipped to window bounds
  awk -v OFS="\t" -v C="$chrom" -v WS="$ws" -v WE="$we" -v WID="$wid" '
    # from regions_with_windows.bed (cols as above)
    $1==C && $5==WID {
      s = ($2>WS ? $2 : WS); e = ($3<WE ? $3 : WE);
      if (e>s) print $1,s,e;
    }' "$tmpdir/regions_with_windows.bed" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$tmpdir/_profile.bed"

  # genomic_bases (from union); fall back to 0 if empty
  GENOMIC_BASES=0
  [[ -s "$tmpdir/_profile.bed" ]] && GENOMIC_BASES=$(awk '{s+=($3-$2)} END{print (s?s:0)}' "$tmpdir/_profile.bed")

  # aligned_bases
  if [[ -s "$tmpdir/_profile.bed" ]]; then
    ALIGNED_BASES=$(
      samtools depth -@ "$cpus" -a -b "$tmpdir/_profile.bed" "$bam" \
      | awk '{sum+=$3} END{print (sum?sum:0)}'
    )
    # I/D/S from CIGAR
    read I_SUM D_SUM S_SUM < <(
      samtools view -@ "$cpus" -F "$SAMTOOLS_FFLAG" -L "$tmpdir/_profile.bed" "$bam" \
      | awk '{
          cig=$6;
          while (match(cig, /([0-9]+)([MIDNSHP=X])/, m)) {
            len=m[1]; op=m[2];
            if (op=="I") I+=len;
            else if (op=="D") D+=len;
            else if (op=="S") S+=len;
            cig=substr(cig, RSTART+RLENGTH);
          }
        } END { printf("%d %d %d\n", (I?I:0), (D?D:0), (S?S:0)) }'
    )
  else
    ALIGNED_BASES=0; I_SUM=0; D_SUM=0; S_SUM=0
  fi
  IDS_SUM=$(( I_SUM + D_SUM + S_SUM ))

  # prefer exact genomic_bases per loop; but if you want the faster bedtools-coverage figure:
  #   GENOMIC_BASES=$(awk -v id="$wid" 'id==$1{print $2; found=1} END{if(!found) print 0}' "$tmpdir/win_genomic_bases.tsv")
  echo -e "${wid}\t${GENOMIC_BASES}\t${ALIGNED_BASES}\t${I_SUM}\t${D_SUM}\t${S_SUM}\t${IDS_SUM}" >> "$tmpdir/window_cov.tsv"
done < "$tmpdir/windows.list"

# Attach BAM covariates to window summary
awk -F'\t' 'NR==FNR{cov[$1]=$0; next}
NR==1{ print $0"\tgenomic_bases\taligned_bases\tI_sum\tD_sum\tS_sum\tIDS_sum"; next }
{
  key=$3
  if (key in cov) {
    split(cov[key], c, "\t")
    print $0"\t"c[2]"\t"c[3]"\t"c[4]"\t"c[5]"\t"c[6]"\t"c[7]
  } else {
    print $0"\t0\t0\t0\t0\t0\t0"
  }
}' "$tmpdir/window_cov.tsv" "$tmpdir/window_summary.tsv" > {3}.${6}.profile.tsv
