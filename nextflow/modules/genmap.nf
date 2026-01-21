process GENMAP {
    def process_name = "genmap"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GenMap/1.3.0-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(ref_genome), path(genome_index_files)
    val(genmap_kmer_length)
    val(genmap_error_tol)
    val(genmap_thresh)

    output: 
    path("genmap_mask.bed"),                                              emit: mask_bed

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${genmap_kmer_length} \
        ${genmap_error_tol} \
        ${genmap_thresh}

    """
}