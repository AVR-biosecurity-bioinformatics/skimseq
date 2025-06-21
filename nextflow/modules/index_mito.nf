process INDEX_MITO {
    def process_name = "index_mito"    
    // tag "-"
    // label "small"
    time '10.m'
    memory '1.GB'
    cpus 1
    // publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "bwa-mem2/2.2.1-GCC-13.3.0:seqtk/1.4-GCC-13.3.0"

    input:
    path(ref_genome)
    val(mito_contig)

    output: 
    tuple path(mito_genome), path("*.fa.*"),        emit: fasta_indexed
    tuple path(mito_bed), path("*.bed"),            emit: mito_bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${ref_genome} \
        ${mito_contig}

    """
}