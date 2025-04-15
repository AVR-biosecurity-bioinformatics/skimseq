process FASTP {
    def process_name = "fastp"    
    // tag "-"
    // label "small"
    time '15.m'
    memory '4.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "fastp/0.23.4-GCC-13.3.0"

    input:
    path(fastq_file)

    output: 
    path("*.trimmed.fastq"),             emit: fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${fastq_file} 

    """
}