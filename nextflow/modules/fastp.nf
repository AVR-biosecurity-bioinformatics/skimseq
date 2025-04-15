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
    tuple val(sample), path(fastq1), path(fastq2)

    output: 
    path("trimmed.R{1,2}.fastq"),             emit: fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${fastq1} \
        ${fastq2}

    """
}