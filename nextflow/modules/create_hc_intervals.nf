process CREATE_HC_INTERVALS {
    def process_name = "create_hc_intervals"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(exclude_bed)
    val(interval_size)
    path(counts)

    output: 
    path("_*.bed"),              emit: interval_bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_size} \
        ${include_bed} \
        "${exclude_bed}" \
        ${ref_fasta} \
        "${counts_files}"

    """
  
}