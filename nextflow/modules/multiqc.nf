process MULTIQC {
    publishDir "${launchDir}/output/modules/multiqc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'

    input:
    path(multiqc_files)
    path(multiqc_config)
    //path(renaming_csv)

    output: 
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots
    
    script:
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash multiqc.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${multiqc_config} 

    """
}