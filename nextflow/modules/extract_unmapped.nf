process EXTRACT_UNMAPPED {
    def process_name = "extract_unmapped"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam), path(bam_index)

    output: 
    tuple val(sample), path(bam), path(bam_index),          emit: bam
    tuple val(sample), path("*.unmapped.R{1,2}.fastq.gz"),  emit: unmapped_fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${bam}

    """
}