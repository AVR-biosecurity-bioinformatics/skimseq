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
