process SUMMARISE_MASKS {
    tag "${ref_genome}"
    publishDir "${launchDir}/output/modules/summarise_masks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'

    input:
    tuple path(ref_genome), path(indexes)
    path(include_bed)
    path(exclude_bed)

    output: 
    path("mask_summary.bed"),              emit: interval_bed
    path("mask_summary.txt"),              emit: summary_file
    
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
        > mask_summary.bed

    # Summarise sequence length by annotation.
    awk 'BEGIN {OFS="\\t"}
        NF >= 4 {sum[\$4] += \$3 - \$2}
        END {
            for (annotation in sum)
                print annotation, sum[annotation]
        }
    ' mask_summary.bed \
        | sort -k1,1 \
        > mask_summary.txt
    """
  
}