process SUMMARISE_MASKS {
    tag "${ref_genome}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(ref_genome), path(indexes)
    path(include_bed)
    path(exclude_bed)

    output: 
    tuple path("mask_summary.bed.gz"), path("mask_summary.bed.gz.tbi"),    emit: mask_summary_bed
    tuple path("mask_pass.bed.gz"), path("mask_pass.bed.gz.tbi"),    emit: mask_pass_bed
    path("mask_summary.txt"),              emit: mask_summary

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    FAI="${ref_genome}.fai"

    # Normalise inputs.
    awk 'BEGIN {OFS="\\t"} NF >= 3 {print \$1, \$2, \$3}' \
        "${include_bed}" \
        | bedtools sort -i - -g "\$FAI" \
        > included_intervals.bed

    awk 'BEGIN {OFS="\\t"} NF >= 3 {
        print \$1, \$2, \$3, (NF >= 4 && \$4 != "" ? \$4 : "Excluded")
    }' \
        "${exclude_bed}" \
        | bedtools sort -i - -g "\$FAI" \
        > excluded_intervals.bed

    # Retain included regions that do not overlap an explicit mask.
    bedtools subtract \
        -a included_intervals.bed \
        -b excluded_intervals.bed \
        | awk 'BEGIN {OFS="\\t"} {print \$1, \$2, \$3, "Included"}' \
        > retained_intervals.bed

    # Label everything outside the explicit masks and retained intervals as excluded.
    cat excluded_intervals.bed retained_intervals.bed \
        | cut -f1-3 \
        | bedtools sort -i - -g "\$FAI" \
        | bedtools merge -i - \
        | bedtools complement -i - -g "\$FAI" \
        | awk 'BEGIN {OFS="\\t"} {print \$1, \$2, \$3, "Excluded"}' \
        > complement_intervals.bed

    # Produce the complete interval classification.
    cat excluded_intervals.bed retained_intervals.bed complement_intervals.bed \
        | bedtools sort -i - -g "\$FAI" \
        | bgzip > mask_summary.bed.gz
    
    tabix -p bed mask_summary.bed.gz

    # Produce a bed of just passing intervals
    zcat mask_summary.bed.gz \
        | awk '\$4 == "Included" {print \$1, \$2, \$3}' OFS='\\t' \
        | bgzip > mask_pass.bed.gz

    tabix -p bed mask_pass.bed.gz

    # Summarise sequence length by annotation.
    zcat mask_summary.bed.gz \
        | awk 'BEGIN {OFS="\\t"}
            NF >= 4 {sum[\$4] += \$3 - \$2}
            END {
                for (annotation in sum)
                    print annotation, sum[annotation]
            }
        ' \
        | sort -k1,1 \
        > mask_summary.txt
    """
  
}