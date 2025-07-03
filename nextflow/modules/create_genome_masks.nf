process CREATE_GENOME_MASKS {
    def process_name = "create_genome_masks"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21:SAMtools/1.21-GCC-13.3.0"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(exclude_bed)
    val(exclude_padding)
    val(mito_contig)
    val(use_reference_hardmasks)
    val(use_reference_softmasks)

    output: 
    path("hard_masked.bed"),              emit: hard_mask
    path("soft_masked.bed"),              emit: soft_mask
    path("mask_summary.bed"),             emit: mask_summary

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
        ${mito_contig} \
        ${ref_fasta} \
        ${use_reference_hardmasks} \
        ${use_reference_softmasks}
        
    """
  
}