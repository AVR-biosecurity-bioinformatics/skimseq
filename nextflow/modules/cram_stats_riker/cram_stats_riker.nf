process CRAM_STATS_RIKER {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/cram_stats_riker", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.txt"),           emit: stats
    tuple val(sample), path("*.pdf"),           emit: plots
    
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