process GENMAP {
    publishDir "${launchDir}/output/modules/genmap", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple path(ref_genome), path(genome_index_files)
    val(genmap_kmer_length)
    val(genmap_error_tol)
    val(genmap_thresh)

    output: 
    path("genmap_mask.bed"),                                              emit: mask_bed

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash genmap.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${genmap_kmer_length} \
        ${genmap_error_tol} \
        ${genmap_thresh}

    """
}