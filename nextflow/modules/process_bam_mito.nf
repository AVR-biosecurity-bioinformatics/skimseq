process PROCESS_BAM_MITO {
    def process_name = "process_bam_mito"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(temp_bam, name: 'temp*.bam')

    output: 
    tuple val(sample), path("sorted.bam"), path("sorted.bam.bai"),        emit: bam
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        "${temp_bam}"

    """
}