process CRAM_STATS {
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/alignment_stats", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.stats.txt"),               emit: stats
    tuple val(sample), path("*.flagstats.txt"),           emit: flagstats
    tuple val(sample), path("*.coverage.txt"),            emit: coverage
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash cram_stats.sh \
        ${task.cpus} \
        ${sample} \
        "${cram}" \
        ${ref_genome}

    """
}