#!/bin/bash
set -euo pipefail

## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = one or more *.stderr.log files

# Collect files exactly as given after $2
LOGS=( "${@:3}" )

PROG_OUT="progress_summary.tsv"
ALLELE_OUT="too_many_alleles.tsv"
SKIP_OUT="skipped_sites.tsv"

# headers
printf "source\tchrom\tpos\telapsed_min\tvariants\tvariants_per_min\n" > "$PROG_OUT"
printf "source\tchrom\tpos\n" > "$ALLELE_OUT"
printf "source\tchrom\tpos\n" > "$SKIP_OUT"

for log in "${LOGS[@]}"; do
  # Source label: basename, drop leading underscore (optional), drop suffix
  base="$(basename -- "$log")"
  src="${base#_}"              # remove leading '_' if present; delete this line to keep it
  src="${src%.stderr.log}"

  # Progress lines
  grep -F 'ProgressMeter -' "$log" 2>/dev/null | \
  awk -v src="$src" '
    /ProgressMeter -/ {
      chrompos=$3; elapsed=$4; nvar=$5; rate=$6
      split(chrompos,a,":"); chrom=a[1]; pos=a[2]
      if (elapsed ~ /^[0-9.]+$/ && nvar ~ /^[0-9]+$/ && rate ~ /^[0-9.]+$/)
        printf("%s\t%s\t%s\t%.3f\t%d\t%.1f\n", src, chrom, pos, elapsed, nvar, rate)
    }' >> "$PROG_OUT" || true

  # Too-many-alleles positions
  grep -F 'has too many alleles in the combined VCF record' "$log" 2>/dev/null | \
  awk -v src="$src" '{
    chrom=""; pos=""
    for(i=1;i<=NF;i++){
      if($i=="Chromosome"){chrom=$(i+1)}
      if($i=="position"){pos=$(i+1)}
    }
    sub(/:$/,"",chrom); gsub(/[^0-9]/,"",pos)
    if(chrom!="" && pos!="") printf("%s\t%s\t%s\n", src, chrom, pos)
  }' >> "$ALLELE_OUT" || true

  # Skipped sites (insufficient data)
  grep -F 'MinimalGenotypingEngine - Some genotypes contained insufficient data' "$log" 2>/dev/null | \
  awk -v src="$src" -F'location ' '
    NF>1 {
      split($2,a,":"); chrom=a[1]; pos=a[2]; gsub(/[[:space:]]/,"",pos)
      if(chrom!="" && pos!="") printf("%s\t%s\t%s\n", src, chrom, pos)
    }' >> "$SKIP_OUT" || true
done

printf "Wrote: %s, %s, %s\n" "$PROG_OUT" "$ALLELE_OUT" "$SKIP_OUT"


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


#if [[ "${DO_PROFILE}" == "true" ]]; then
#    # Wall-clock timer end
#    END_NS=$(date +%s%N)
#    END_TS="$(date -u +%FT%TZ)"
#    ELAPSED_SEC=$(( (END_NS - START_NS)/1000000000 ))

    # Output TSV header 
#    PROFILE_TSV="${3}.${6}.profile.tsv"
#    echo -e "sample\tinterval_hash\tintervals\tgenomic_bases\taligned_bases\ti_sum\td_sum\ts_sum\tids_sum\tstart_iso8601\tend_iso8601\telapsed_seconds\thc_threads\thc_min_pruning\thc_min_dangling\thc_max_reads_startpos\thc_rmdup\thc_minbq\thc_minmq\tploidy" > "${PROFILE_TSV}"

    # Append profiling row
#    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
#      "${3}" "${6}" "${INTERVALS_COUNT}" "${GENOMIC_BASES}" "${ALIGNED_BASES}" \
#      "${I_SUM}" "${D_SUM}" "${S_SUM}" "${IDS_SUM}" \
#      "${START_TS}" "${END_TS}" "${ELAPSED_SEC}" \
#      "${1}" "${11}" "${12}" "${13}" "${14}" "${15}" "${16}" "${17}" \
#      >> "${PROFILE_TSV}"
#fi
