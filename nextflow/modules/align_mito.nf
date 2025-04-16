process ALIGN_MITO {
    def process_name = "align_mito"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BWA/0.7.18-GCCcore-13.3.0:SAMtools/1.21-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2), path(json)
    tuple path(mito_genome), path(mito_index_files)

    output: 
    tuple val(sample), path("*.tempmito.bam"),             emit: bam
    
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
        ${json} \
        ${mito_genome}

    """
}