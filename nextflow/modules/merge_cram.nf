process MERGE_CRAM {
    publishDir "${launchDir}/output/modules/merge_cram", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir params.cram_store, mode: 'copy', pattern: "*.cram*", enabled: "${ params.output_cram ? true : false }"

    input:
    tuple val(sample), val(lib), path(cram) 
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram
    tuple val(sample), path("*.markdup.json"),                                emit: markdup

    script:
    """
    #!/usr/bin/env bash

    # Write list of cram files to process
    printf "%s\n" ${cram} > cram.list
    
    ### run process script
    bash merge_cram.sh \
        ${task.cpus} \
        ${sample} \
        ${ref_genome}

    """
}