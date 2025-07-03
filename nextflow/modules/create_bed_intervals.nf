process CREATE_BED_INTERVALS {
    def process_name = "create_bed_intervals"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(hard_masks_bed)
    path(soft_masks_bed)
    val(interval_n)
    val(interval_size)
    val(interval_subdivide)
    val(interval_include_hard_masks)
    val(interval_include_soft_masks)
    
    output: 
    path("*.bed"),              emit: interval_list
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${interval_n} \
        ${interval_size} \
        ${include_bed} \
        ${hard_masks_bed} \
        ${soft_masks_bed} \
        ${interval_include_hard_masks} \
        ${interval_include_soft_masks} \
        ${interval_subdivide} \
        ${ref_fasta}

        
    """
  
}