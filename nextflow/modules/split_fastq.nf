process SPLIT_FASTQ {
    def process_name = "split_fastq"    
    // tag "-"
    // label "small"
    time '10.m'
    memory '4.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "seqtk/1.4-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2), path(json)
    val(chunk_size)

    output: 
    tuple val(sample), path("R1.*.fastq.gz"), path("R2.*.fastq.gz"), path(json),     emit: fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${fastq1} \
        ${fastq2} \
        ${chunk_size}

    """
}