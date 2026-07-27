process EXTRACT_GENOME_MASKS {
    tag "${ref_genome}"
    publishDir "${launchDir}/output/modules/extract_genome_masks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple path(ref_genome), path(indexes)
    path(include_bed)
    path(exclude_bed)
    val(exclude_padding)
    val(use_reference_hardmasks)
    val(use_reference_softmasks)

    output: 
    path("genome_masks.bed"),                    emit: mask_bed

    script:
    """
    #!/usr/bin/env bash
    
    # Included intervals is just the reference genome bed if specific intervals were not provided
    # Just keep the important columns
    cat ${include_bed} | cut -f1-3 > included_intervals.bed

    # Ensure the intermediate and final mask files exist even when no masks are enabled or no mask intervals are found.
    : > concat_masks.bed
    : > genome_masks.bed

    # Add any excluded intervals to genome mask bed if file is not empty
    if [ -s ${exclude_bed} ] ; then  
        # Keep just the important columns
        awk ' BEGIN { OFS = "\\t" } NF >= 3 { print \$1, \$2, \$3, "Excluded" } ' "${exclude_bed}" \
        >> concat_masks.bed
    fi

    # If use_reference_hardmasks is true, add uppercase N regions from the reference.
    if [[ "${use_reference_hardmasks}" == "true" ]]; then
        seqkit locate \
            --only-positive-strand \
            --use-regexp \
            --pattern "N+" \
            --bed \
            --id-regexp "^(\\\\S+)" \
            "${ref_genome}" \
            | awk ' BEGIN { OFS = "\\t" } NF >= 3 { print \$1, \$2, \$3, "NRef" } ' \
            >> concat_masks.bed
    fi

    # If exclude_reference_genome_softmasks is true, add lowercase regions from the reference.
        if [ ${use_reference_softmasks} == "true" ] ; then
            seqkit locate \
                --only-positive-strand \
                --use-regexp \
                --pattern "[a-z]+" \
                --bed \
                --id-regexp "^(\\\\S+)" \
                "${ref_genome}" \
                | awk 'BEGIN { OFS = "\\t" } NF >= 3 { print \$1, \$2, \$3, "SoftMaskRef" }' \
                >> concat_masks.bed
    fi

    # Sort and merge all overlapping or book-ended mask intervals.
    if [[ -s concat_masks.bed ]]; then
        bedtools sort \
            -i concat_masks.bed |
            bedtools merge \
                -i - \
                -c 4 \
                -o distinct \
            > genome_masks.bed
    fi    
    """
  
}