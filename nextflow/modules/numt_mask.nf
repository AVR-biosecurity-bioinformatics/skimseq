process NUMT_MASK {
    publishDir "${launchDir}/output/modules/numt_mask", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple path(ref_genome), path(genome_index_files)
    tuple path(mito_genome), path(mito_index_files)
    val(numt_min_length)
    val(numt_max_gap)

    output: 
    path("numt_mask.bed"),                                              emit: mask_bed

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash numt_mask.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${mito_genome} \
        ${numt_min_length} \
        ${numt_max_gap}

    """
}