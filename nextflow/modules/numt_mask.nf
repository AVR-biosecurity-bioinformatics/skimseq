process NUMT_MASK {
    def process_name = "numt_mask"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "MUMmer/4.0.1-GCCcore-13.3.0:SeqKit/2.8.2:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(ref_genome), path(genome_index_files)
    tuple path(mito_genome), path(mito_index_files)
    val(numt_min_length)
    val(numt_max_gap)

    output: 
    path("numt_mask.bed"),                                              emit: mask_bed

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${mito_genome} \
        ${numt_min_length} \
        ${numt_max_gap}

    """
}