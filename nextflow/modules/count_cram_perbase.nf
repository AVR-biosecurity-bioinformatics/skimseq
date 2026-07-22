process COUNT_CRAM_PERBASE {
    // tag "-"
    publishDir "${launchDir}/output/modules/count_cram_perbase", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/alignment_stats", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    val(hc_rmdup)
    val(hc_minbq)
    val(hc_minmq)
    
    output: 
    tuple val(sample), path("${sample}.perbase.bed.gz"),  path("${sample}.perbase.bed.gz.tbi"),   emit: perbase

    script:
    """
    #!/usr/bin/env bash

    ### run process script
    bash count_cram_perbase.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${cram}" \
        ${sample} \
        ${ref_genome} \
        ${hc_rmdup} \
        ${hc_minbq} \
        ${hc_minmq}
        
    """
}
