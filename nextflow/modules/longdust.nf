process LONGDUST {
    publishDir "${launchDir}/output/modules/longdust", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple path(ref_genome), path(genome_index_files)
    val(longdust_kmer_length)
    val(longdust_window_size)
    val(longdust_thresh)

    output: 
    path("longdust_mask.bed"),                                              emit: mask_bed

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash longdust.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${longdust_kmer_length} \
        ${longdust_window_size} \
        ${longdust_thresh}

    """
}