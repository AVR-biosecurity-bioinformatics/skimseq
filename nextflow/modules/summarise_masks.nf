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
    
    # Normalise and sort the included analysis intervals.
    awk 'BEGIN { OFS = "\\t" } NF >= 3 { print \$1, \$2, \$3 }' "${include_bed}" \
        | bedtools sort -i - -g "${ref_genome}.fai" \
        > included_intervals.bed

    # First add the masked intervals
    cp ${exclude_bed} all_intervals.bed

    # Then add any included intervals that arent in the masks
    bedtools sort -i included_intervals.bed -g ${ref_genome}.fai \
        | bedtools subtract -a stdin -b ${exclude_bed} \
        | cut -f1-4 >> all_intervals.bed

    # Then add any remaining intervals of the genome that werent in the exclude or include masks
    bedtools sort -i all_intervals.bed -g ${ref_genome}.fai \
        | bedtools complement -i stdin -g ${ref_genome}.fai \
        | sed 's/\s*$/\tExcluded/' >> all_intervals.bed

    # Create sorted mask summary file
    bedtools sort -i all_intervals.bed -g ${ref_genome}.fai > mask_summary.bed

    # Tabulate by annotation type
    awk '{len=\$3-\$2; sum[\$4]=sum[\$4]+len} END {for (anno in sum) print anno, sum[anno]}'  mask_summary.bed > mask_summary.txt

    """
  
}