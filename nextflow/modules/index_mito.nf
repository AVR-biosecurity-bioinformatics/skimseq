process INDEX_MITO {
    publishDir "${launchDir}/output/modules/index_mito", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ref_genome)
    val(mito_contig)

    output: 
    tuple path("*.fa"), path("*.{fa.*,fna.*,dict}"),    emit: fasta_indexed
    path("mito.bed"),                                   emit: bed
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash index_mito.sh \
        ${task.cpus} \
        ${ref_genome} \
        "${mito_contig}"

    """
}