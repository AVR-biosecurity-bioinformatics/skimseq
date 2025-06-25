process CREATE_BED_INTERVALS {
    def process_name = "create_bed_intervals"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    tuple path(ref_fasta), path(indexes)
    val(interval_size)
    path(interval_bed)
    val(interval_padding)
    path(exclude_bed)
    val(exclude_padding)
    val(mito_contig)

    output: 
    path("*.bed"),              emit: interval_list
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${interval_size} \
        ${included_intervals} \
        ${interval_padding} \
        ${exclude_intervals} \
        ${exclude_padding} \
        ${mito_contig}

    """
  
}