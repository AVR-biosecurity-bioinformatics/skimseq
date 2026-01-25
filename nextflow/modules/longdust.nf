process LONGDUST {
    def process_name = "longdust"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "longdust/1.4-GCCcore-13.3.0:SeqKit/2.8.2:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(ref_genome), path(genome_index_files)
    val(longdust_kmer_length)
    val(longdust_window_size)
    val(longdust_thresh)

    output: 
    path("longdust_mask.bed"),                                              emit: mask_bed

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${longdust_kmer_length} \
        ${longdust_window_size} \
        ${longdust_thresh}

    """
}