process PROFILE_HC {
    def process_name = "profile_hc"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(sample), val(interval_hash), path(logfile), path(assembly_regions)

    output: 
    path("*.progress_summary.tsv"),                                             emit: summary

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        ${interval_hash} \
        ${logfile}

    """
}