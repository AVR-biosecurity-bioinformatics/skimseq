process MULTIQC {
    def process_name = "multiqc"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "MultiQC/1.28-foss-2024a"

    input:
    path(multiqc_files)

    output: 
    path "*multiqc_report.html", emit: report
    //path "*_data"              , emit: data
    //path "*_plots"             , emit: plots
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} 

    """
}