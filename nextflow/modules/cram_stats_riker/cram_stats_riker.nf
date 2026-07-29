process CRAM_STATS {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
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
    set -euo pipefail

    riker multi \
        --threads ${task.cpus} \
        -i "${cram}" \
        -r ${ref_genome} \
        -o ${sample} \
        --tools alignment isize basic gcbias wgs error \
        --error::stratify-by read_num,cycle bq


    """
}