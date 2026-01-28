process CREATE_INTERVAL_CHUNKS {
    def process_name = "create_interval_chunks"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BEDOPS/2.4.42-foss-2024a:parallel/20240722-GCCcore-13.3.0:BCFtools/1.22-GCC-13.3.0"

    input:
    tuple val(sample), path(counts_bed), path(counts_tbi)
    tuple path(ref_genome), path(genome_index_files)
    path(contig_bed)
    val(counts_per_chunk)
    val(min_interval_gap)
    val(split_large_intervals)
    val(include_zero)

    output: 
    tuple val(sample), path("_*.bed"),              emit: interval_bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Write list of bed files to process
    printf "%s\n" ${counts_bed} > counts_files.list
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${counts_per_chunk} \
        ${split_large_intervals} \
        ${min_interval_gap} \
        ${contig_bed} \
        ${include_zero} 

    """
  
}