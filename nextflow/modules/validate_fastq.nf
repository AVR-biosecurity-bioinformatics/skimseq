process VALIDATE_FASTQ {
    def process_name = "validate_fastq"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "seqtk/1.4-GCC-13.3.0"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), path(fastq1), path(fastq2), stdout, emit: fastq_with_status
    
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