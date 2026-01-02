process PROFILE_HC {
    def process_name = "profile_hc"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0"

    input:
    tuple val(sample), val(interval_hash), path(cram), path(cram_index), path(gvcf), path(gvcf_index), path(logfile), path(assembly_regions)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    path("*.profile.tsv"),                                             emit: summary

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        ${interval_hash} \
        ${logfile} \
        ${assembly_regions} \
        ${cram} \
        ${gvcf} \
        ${ref_genome}

    """
}