process CREATE_BED_INTERVALS {
    def process_name = "create_bed_intervals"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21"

    input:
    tuple path(ref_fasta), path(indexes)
    val(interval_splitn)
    val(interval_n)
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
        ${interval_splitn} \
        ${interval_n} \
        ${interval_bed} \
        ${interval_padding} \
        ${exclude_bed} \
        ${exclude_padding} \
        ${mito_contig} \
        ${ref_fasta}
    """
  
}