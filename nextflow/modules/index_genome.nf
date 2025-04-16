process INDEX_GENOME {
    def process_name = "index_genome"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    // publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BWA/0.7.18-GCCcore-13.3.0"

    input:
    path(ref_genome)

    output: 
    tuple path(ref_genome), path("*.fa.*"),             emit: fasta_indexed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${ref_genome}

    """
}