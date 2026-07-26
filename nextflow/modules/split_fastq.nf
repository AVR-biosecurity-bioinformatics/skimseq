process SPLIT_FASTQ {
    tag "${lib}"
    publishDir "${launchDir}/output/modules/split_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)
    val(chunk_size)

    output: 
    tuple val(sample), val(lib), path("intervals_${lib}.csv"), emit: fastq_interval 
    tuple val(sample), val(lib), path("nchunks_${lib}.txt"), emit: nchunks

    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash split_fastq.sh \
        ${task.cpus} \
        ${lib} \
        ${fastq1} \
        ${fastq2} \
        ${chunk_size}

    """
}