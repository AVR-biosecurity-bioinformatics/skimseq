process COUNT_VCF_BED {
    def process_name = "count_vcf_bed"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(sample), path(gvcf), path(tbi)
    path(interval_bed)
    path(exclude_bed)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*counts.bed"),                 emit: counts
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${gvcf}" \
        ${ref_genome} \
        ${interval_bed} \
        ${exclude_bed} \
        ${sample}
        
    """
}
