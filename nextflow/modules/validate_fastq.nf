process VALIDATE_FASTQ {
    def process_name = "validate_fastq"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SeqKit/2.8.2"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), stdout, emit: status
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${lib} \
        ${fastq1} \
        ${fastq2}

    """
}