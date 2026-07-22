process SUM_COVERED_INTERVALS {
    def process_name = "sum_covered_intervals"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(count_bed),  path(tbi)
    path(exclude_bed)
    
    output: 
    tuple val(sample), path("${sample}.covered.bed.gz"),  path("${sample}.covered.bed.gz.tbi"),   emit: counts

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        ${count_bed} \
        ${exclude_bed}
    """
}
