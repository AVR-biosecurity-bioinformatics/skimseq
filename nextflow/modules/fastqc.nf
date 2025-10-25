process FASTQC {
    def process_name = "fastqc"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "FastQC/0.12.1-Java-11"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)

    output: 
    path("*.zip"),             emit: results
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${fastq1} \
        ${fastq2} \
        ${sample} 

    """
}