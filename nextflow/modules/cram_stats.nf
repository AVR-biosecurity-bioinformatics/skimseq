process CRAM_STATS {
    publishDir "${launchDir}/output/modules/cram_stats", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
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

    # Output sample coverage statistics
    samtools coverage \
        -Q1 \
        -q1 \
        --reference ${ref_genome} \
        "${cram}" \
        > ${sample}.coverage.txt

    # Output flag statistics
    samtools flagstats \
        --threads ${task.cpus} \
        "${cram}" \
        > ${sample}.flagstats.txt

    # Output comprehensive statistics
    samtools stats \
        --threads ${task.cpus} \
        --reference ${ref_genome} \
        "${cram}" \
        > ${sample}.stats.txt

    """
}