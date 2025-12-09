process CREATE_INTERVAL_CHUNKS_HC {
    def process_name = "create_interval_chunks_hc"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(sample), path(counts_files)
    val(counts_per_chunk)
    val(split_overweight)

    output: 
    tuple val(sample), path("_*.bed"),              emit: interval_bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${counts_per_chunk} \
        ${split_overweight} \
        "${counts_files}" \


    """
  
}