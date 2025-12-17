process CREATE_INTERVAL_CHUNKS_JC {
    def process_name = "create_interval_chunks_jc"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:parallel/20240722-GCCcore-13.3.0:BCFtools/1.22-GCC-13.3.0:"

    input:
    tuple path(counts_bed), path(counts_tbi)
    val(counts_per_chunk)
    val(split_large_intervals)
    val(min_interval_gap)
    tuple path(ref_genome), path(genome_index_files)

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
        ${ref_genome} \
        ${counts_per_chunk} \
        ${split_large_intervals} \
        ${min_interval_gap} \
        "${counts_bed}" 

    """
  
}