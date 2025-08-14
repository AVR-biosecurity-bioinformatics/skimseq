process INDEX_MITO {
    def process_name = "index_mito"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "bwa-mem2/2.2.1-GCC-13.3.0:SAMtools/1.21-GCC-13.3.0:seqtk/1.4-GCC-13.3.0"

    input:
    path(ref_genome)
    val(mito_contig)

    output: 
    tuple path("*.fa"), path("*.{fa.*,fna.*,dict}"),    emit: fasta_indexed
    path("*.bed"),                                      emit: bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${ref_genome} \
        "${mito_contig}"

    """
}