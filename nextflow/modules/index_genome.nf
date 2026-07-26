process INDEX_GENOME {
    tag "${ref_genome}"
    publishDir "${launchDir}/output/modules/index_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ref_genome)
    val(min_chr_length)
    
    output: 
    tuple path(ref_genome), path("*.{fa.*,fna.*,fasta.*,dict}"),     emit: fasta_indexed
    path("genome.bed"),                                              emit: genome_bed
    path("long.bed"),                                                emit: long_bed
    path("short.bed"),                                               emit: short_bed

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash index_genome.sh \
        ${task.cpus} \
        ${ref_genome} \
        ${min_chr_length}

    """
}