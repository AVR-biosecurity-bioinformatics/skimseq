process MULTIQC {
    def process_name = "multiqc"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "MultiQC/1.28-foss-2024a"

    input:
    path(multiqc_files)
    path(renaming_csv)
    path(multiqc_config)

    output: 
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${multiqc_config} \
        ${renaming_csv}

    """
}