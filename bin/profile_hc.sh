#!/bin/bash
set -euo pipefail

## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample
# $4 = interval_hash
# $5 = .stderr.log file

# Header with the two new columns first
printf "sample\tinterval_hash\tchrom\tpos\telapsed_min\tregions\tregions_per_min\n" > progress_summary.tsv

# Stream: grep only ProgressMeter lines, awk extracts data rows and prefixes sample/hash
grep -F 'ProgressMeter -' "${5}" \
    | awk -v s="${3}" -v h="${4}" '
    # match data rows like:
    # ... ProgressMeter -   CM028320.1:49552   0.0   24   534.1
    match($0, /ProgressMeter -[[:space:]]*([[:alnum:]_.]+):([0-9]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9.]+)/, m) {
        printf("%s\t%s\t%s\t%s\t%.3f\t%d\t%.1f\n", s, h, m[1], m[2], m[3], m[4], m[5])
    } ' >> progress_summary.tsv


# Transform assembly regions output into a bed file
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 || $1 ~ /^#/ { next }             # skip track/header
    $4 == "end-marker" { next }             # drop 0-length markers
    $4 ~ /^size=/ {
    size = $4; sub(/^size=/,"",size)
    v = $5 + 0
    active = (v > 0 ? 1 : 0)              # +1 active, -1 inactive in GATK IGV track
    print $1, $2, $3, active, size
    }' activity_regions.tsv \
    | sort -k1,1 -k2,2n > hc_ar.bed


# Turn ticks into spans between consecutive ticks
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 { next }  # skip header
{
  s=$1; h=$2; c=$3; p=$4+0; em=$5+0; r=$6+0; rpm=$7+0
  key=s"|"h"|"c
  if(!(key in seen)) {
    # seed a tiny 1bp stub for the first tick per (sample,hash,chrom)
    start0=(p>0?p-1:0); end=p; if(end<=start0){end=start0+1}
    seglen=end-start0
    d_em=em; d_r=r
    print c,start0,end,s,h,p,em,d_em,r,d_r,rpm,seglen,(d_em*60.0)/(seglen/1000.0),(d_r)/(seglen/1000.0)
    seen[key]=1; P[key]=p; EM[key]=em; R[key]=r
    next
  }
  start0=(P[key]>0?P[key]-1:0); end=p; if(end<=start0){end=start0+1}
  seglen=end-start0
  d_em=em-EM[key]; if(d_em<0)d_em=0
  d_r=r-R[key];   if(d_r<0)d_r=0
  print c,start0,end,s,h,p,em,d_em,r,d_r,rpm,seglen,(d_em*60.0)/(seglen/1000.0),(d_r)/(seglen/1000.0)
  P[key]=p; EM[key]=em; R[key]=r
}' progress_summary.tsv \
| awk 'BEGIN{OFS="\t"} $3>$2' \
| sort -k1,1 -k2,2n \
> hc_ticks.bed


# Attribute runtime to assembly regions
# Intersect; -wo appends overlap length as last field
bedtools intersect -a hc_ar.bed -b hc_ticks.bed -wo \
> ar_x_ticks.tmp

# Scale each tick delta to the overlapped fraction of the tick span
awk -F'\t' -v OFS='\t' '
{
  # a: hc_ar.bed -> $1..$5   (chrom,start,end,active,size)
  # b: hc_ticks.bed -> $6..$19
  chrom=$1; astart=$2; aend=$3; active=$4; asize=$5
  bchrom=$6; bstart=$7; bend=$8
  sample=$9; hash=$10
  end_pos=$11; end_em=$12; d_em=$13; end_r=$14; d_r=$15; end_rpm=$16
  seg_len=$17; sec_per_kb=$18; reg_per_kb=$19
  ov=$20

  if (ov<=0 || seg_len<=0) next

  frac = ov / seg_len
  d_em_clip = d_em * frac
  d_r_clip  = d_r  * frac

  # Emit per-overlap row (can aggregate later)
  print chrom, astart, aend, active, asize, sample, hash, ov, d_em_clip, d_r_clip
}' ar_x_ticks.tmp \
> ar_alloc.tsv

# Calculate covariates for assembly regions


# # ------------- PROFILING CODE  -------------
#    # Build bed for profiling if provided
#    if [[ -n "${9}" && "${9}" != "none" && -s "${9}" ]]; then
#        bedtools slop -i "${9}" -g "${5}.fai" -b "${10}" \
#        | bedtools subtract -a "${7}" -b stdin \
#        | bedtools merge > profile.bed
#    else
#        bedtools merge -i "${7}" > profile.bed
#    fi
#
#    # Count intervals and bases
#    INTERVALS_COUNT=$(wc -l < profile.bed)
#    GENOMIC_BASES=$(awk '{s+=($3-$2)} END{print (s? s:0)}' profile.bed)
#
#    # Duplicate handling mirrors HC NotDuplicateReadFilter switch
#    if [[ ${14} == "false" ]];    then SAMTOOLS_FFLAG=260; else SAMTOOLS_FFLAG=1284; fi
#
#    if [[ "${INTERVALS_COUNT}" -eq 0 || "${GENOMIC_BASES}" -eq 0 ]]; then
#        ALIGNED_BASES=0
#        I_SUM=0; D_SUM=0; S_SUM=0
#    else
#        # Sum of per-base depth across intervals
#        ALIGNED_BASES=$(
#            samtools depth -@ "${1}" -a -b profile.bed "${4}" \
#            | awk '{sum+=$3} END{print (sum? sum:0)}'
#        )
#
#        # Separate I / D / S sums from CIGARs
#        read I_SUM D_SUM S_SUM < <(
#            samtools view -@ "${1}" -F "${SAMTOOLS_FFLAG}" -L profile.bed "${4}" \
#            | awk '{
#                cig=$6; i_sum=0; d_sum=0; s_sum=0
#                while (match(cig, /([0-9]+)([MIDNSHP=X])/, m)) {
#                    len=m[1]; op=m[2]
#                    if (op=="I") i_sum+=len
#                    else if (op=="D") d_sum+=len
#                    else if (op=="S") s_sum+=len
#                   cig=substr(cig, RSTART+RLENGTH)
#                }
#                I+=i_sum; D+=d_sum; S+=s_sum
#            } END { printf("%d %d %d\n", (I?I:0), (D?D:0), (S?S:0)) }'
#        )
#    fi
#
#    IDS_SUM=$(( I_SUM + D_SUM + S_SUM ))
