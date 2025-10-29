process SPLIT_FASTQ {
    def process_name = "split_fastq"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "seqtk/1.4-GCC-13.3.0"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)
    val(chunk_size)

    output: 
    tuple val(sample), val(lib), path("intervals_${lib}.csv"), emit: fastq_interval 

    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${lib} \
        ${fastq1} \
        ${fastq2} \
        ${chunk_size}

    """
}