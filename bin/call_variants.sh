#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample name
# $4 = cram file
# $5 = ref genome
# $6 = interval hash
# $7 = interval_bed
# $8 = interval_padding
# $9 = exclude_bed
# $10 = exclude_padding
# $11 = hc_min_pruning
# $12 = hc_min_dangling_length
# $13 = hc_max_reads_startpos
# $14 = hc_rmdup
# $15 = hc_minbq
# $16 = hc_minmq
# $17 = ploidy
# $18 = profile_gatk

# 1GB of memory should be retained outside the java heap
java_mem=$(( ( ${2} - 1 )))

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# parse filtering options as flags
if [[ ${14} == "false" ]];    then RMDUP="-DF NotDuplicateReadFilter";                  else RMDUP=""; fi

# Create list of crams to be processed
echo ${4} | tr ' ' '\n' > cram.list

 # ------------- PROFILING CODE  -------------

# Toggle runtime profiling from arg $18
DO_PROFILE="${18:-false}"

if [[ "${DO_PROFILE}" == "true" ]]; then
    # Build bed for profiling if provided
    if [[ -n "${9}" && "${9}" != "none" && -s "${9}" ]]; then
        bedtools slop -i "${9}" -g "${5}.fai" -b "${10}" \
        | bedtools subtract -a "${7}" -b stdin \
        | bedtools merge > profile.bed
    else
        bedtools merge -i "${7}" > profile.bed
    fi

    # Count intervals and bases
    INTERVALS_COUNT=$(wc -l < profile.bed)
    GENOMIC_BASES=$(awk '{s+=($3-$2)} END{print (s? s:0)}' profile.bed)

    # Duplicate handling mirrors HC NotDuplicateReadFilter switch
    if [[ ${14} == "false" ]];    then SAMTOOLS_FFLAG=260; else SAMTOOLS_FFLAG=1284; fi

    if [[ "${INTERVALS_COUNT}" -eq 0 || "${GENOMIC_BASES}" -eq 0 ]]; then
        ALIGNED_BASES=0
        I_SUM=0; D_SUM=0; S_SUM=0
    else
        # Sum of per-base depth across intervals
        ALIGNED_BASES=$(
            samtools depth -@ "${1}" -a -b profile.bed "${4}" \
            | awk '{sum+=$3} END{print (sum? sum:0)}'
        )

        # Separate I / D / S sums from CIGARs
        read I_SUM D_SUM S_SUM < <(
            samtools view -@ "${1}" -F "${SAMTOOLS_FFLAG}" -L profile.bed "${4}" \
            | awk '{
                cig=$6; i_sum=0; d_sum=0; s_sum=0
                while (match(cig, /([0-9]+)([MIDNSHP=X])/, m)) {
                    len=m[1]; op=m[2]
                    if (op=="I") i_sum+=len
                    else if (op=="D") d_sum+=len
                    else if (op=="S") s_sum+=len
                    cig=substr(cig, RSTART+RLENGTH)
                }
                I+=i_sum; D+=d_sum; S+=s_sum
            } END { printf("%d %d %d\n", (I?I:0), (D?D:0), (S?S:0)) }'
        )
    fi

    IDS_SUM=$(( I_SUM + D_SUM + S_SUM ))

    # Wall-clock timer start
    START_TS="$(date -u +%FT%TZ)"
    START_NS=$(date +%s%N)
fi

# call variants by sample * interval chunk
# NOTE: need to use assembly region padding rather than interval_padding to avoid overlapping variants
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=${1}" HaplotypeCaller \
    -R $5 \
    -I cram.list \
    -L $7 \
    --native-pair-hmm-threads ${1} \
    --assembly-region-padding ${8} \
    --exclude-intervals ${9} \
    --interval-exclusion-padding ${10} \
    --interval-merging-rule ALL \
    --min-pruning ${11} \
    --min-dangling-branch-length ${12} \
    --max-reads-per-alignment-start ${13} \
    $RMDUP \
    --min-base-quality-score ${15} \
    --minimum-mapping-quality ${16} \
    --mapping-quality-threshold-for-genotyping ${16} \
    -ploidy ${17} \
    -ERC GVCF \
    -O tmp.g.vcf.gz

# NOTE: Haplotypecaller ALWAYS outputs intervals in the GVCF, even if there are no reads - so drop these with bcftools
bcftools view \
    -e 'ALT="<NON_REF>" && (MAX(FORMAT/DP)=0 || MAX(FORMAT/MIN_DP)=0 || MAX(FORMAT/GQ)=0)' \
    -Oz -o ${3}.${6}.g.vcf.gz tmp.g.vcf.gz
bcftools index -t ${3}.${6}.g.vcf.gz

if [[ "${DO_PROFILE}" == "true" ]]; then
    # Wall-clock timer end
    END_NS=$(date +%s%N)
    END_TS="$(date -u +%FT%TZ)"
    ELAPSED_SEC=$(( (END_NS - START_NS)/1000000000 ))

    # Output TSV header 
    PROFILE_TSV="${3}.${6}.profile.tsv"
    echo -e "sample\tinterval_hash\tintervals\tgenomic_bases\taligned_bases\ti_sum\td_sum\ts_sum\tids_sum\tstart_iso8601\tend_iso8601\telapsed_seconds\thc_threads\thc_min_pruning\thc_min_dangling\thc_max_reads_startpos\thc_rmdup\thc_minbq\thc_minmq\tploidy" > "${PROFILE_TSV}"

    # Append profiling row
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${3}" "${6}" "${INTERVALS_COUNT}" "${GENOMIC_BASES}" "${ALIGNED_BASES}" \
      "${I_SUM}" "${D_SUM}" "${S_SUM}" "${IDS_SUM}" \
      "${START_TS}" "${END_TS}" "${ELAPSED_SEC}" \
      "${1}" "${11}" "${12}" "${13}" "${14}" "${15}" "${16}" "${17}" \
      >> "${PROFILE_TSV}"
fi

rm -f tmp*
