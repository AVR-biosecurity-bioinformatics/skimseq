process CREATE_GENOME_MASKS {
    def process_name = "create_genome_masks"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21:SAMtools/1.21-GCC-13.3.0"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(excluded_bed)
    val(excluded_padding)
    val(mitochondrial_contig)
    val(exclude_reference_hardmasks)
    val(exclude_reference_softmasks)
    
    output: 
    path("hard_masked.bed"),              emit: hard_mask
    path("soft_masked.bed"),              emit: soft_mask
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${include_bed} \
        ${excluded_bed} \
        ${excluded_padding} \
        ${mitochondrial_contig} \
        ${ref_fasta} \
        ${exclude_reference_hardmasks} \
        ${exclude_reference_softmasks}
        
    """
  
}