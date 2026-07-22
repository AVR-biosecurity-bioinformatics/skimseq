process EXTRACT_GENOME_MASKS {
    def process_name = "extract_genome_masks"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

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