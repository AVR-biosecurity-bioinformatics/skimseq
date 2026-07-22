process EXTRACT_GENOME_MASKS {
    // tag "-"
    publishDir "${launchDir}/output/modules/extract_genome_masks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

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
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash extract_genome_masks.sh \
        ${task.cpus} \
        ${include_bed} \
        ${exclude_bed} \
        ${exclude_padding} \
        ${use_reference_hardmasks} \
        ${use_reference_softmasks} \
        ${ref_fasta}
        
    """
  
}