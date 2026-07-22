process SUM_COVERED_INTERVALS {
    publishDir "${launchDir}/output/modules/sum_covered_intervals", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(count_bed),  path(tbi)
    path(exclude_bed)
    
    output: 
    tuple val(sample), path("${sample}.covered.bed.gz"),  path("${sample}.covered.bed.gz.tbi"),   emit: counts

    script:
    """
    #!/usr/bin/env bash

    ### run process script
    bash sum_covered_intervals.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        ${count_bed} \
        ${exclude_bed}
    """
}
