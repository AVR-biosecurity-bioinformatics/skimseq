process CONSENSUS_MITO {
    def process_name = "consensus_mito"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam), path(bam_index)
    tuple path(mito_genome), path(mito_index_files)

    output: 
    tuple val(sample), path("*.mito.fa.gz"),        emit: fasta
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${bam} \
        ${mito_genome}

    """
}