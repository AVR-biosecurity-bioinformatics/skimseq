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
    
    ### run process script
    bash summarise_masks.sh \
        ${task.cpus} \
        ${include_bed} \
        ${exclude_bed} \
        ${ref_genome}

    """
  
}