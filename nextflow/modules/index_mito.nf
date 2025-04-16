process INDEX_MITO {
    def process_name = "index_mito"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BWA/0.7.18-GCCcore-13.3.0"

    input:
    path(mito_genome)

    output: 
    tuple path(mito_genome), path("*.fa.*"),             emit: fasta_indexed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${mito_genome}

    """
}