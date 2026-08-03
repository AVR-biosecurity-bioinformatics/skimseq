process CREATE_INTERVAL_CHUNKS {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/create_interval_chunks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(counts_bed), path(counts_tbi)
    tuple path(ref_genome), path(genome_index_files)
    path(include_bed)
    val(counts_per_chunk)
    val(min_interval_gap)
    val(split_large_intervals)
    val(include_zero)

    output: 
    tuple val(sample), path("*.bed.gz"), path("*.bed.gz.tbi"),    emit: interval_bed
    
    script:
    """
    #!/usr/bin/env bash

    # Write list of bed files to process
    printf "%s\n" ${counts_bed} > counts_files.list
    
    ### run process script
    bash create_interval_chunks.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${counts_per_chunk} \
        ${split_large_intervals} \
        ${min_interval_gap} \
        ${include_bed} \
        ${include_zero} 

    """
  
}