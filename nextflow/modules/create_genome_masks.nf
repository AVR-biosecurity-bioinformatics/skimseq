process CREATE_GENOME_MASKS {
    def process_name = "create_genome_masks"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:SeqKit/2.8.2"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(exclude_bed)
    val(exclude_padding)
    val(use_reference_hardmasks)
    val(use_reference_softmasks)

    output: 
    path("genome_masks.bed"),                    emit: mask_bed

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${include_bed} \
        ${exclude_bed} \
        ${exclude_padding} \
        ${use_reference_hardmasks} \
        ${use_reference_softmasks} \
        ${ref_fasta}
        
    """
  
}