process FASTP {
    def process_name = "fastp"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '4.GB'
    cpus 4
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "fastp/0.23.4-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), path("trimmed.R1.fastq"), path("trimmed.R2.fastq"), path("*.json"),     emit: fastq
    
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