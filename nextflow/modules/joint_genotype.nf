process JOINT_GENOTYPE {
    def process_name = "joint_genotype"    
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_bed), path(genomicsdb)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)
    val(exclude_padding)
    val(output_invariant)

    output: 
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${genomicsdb} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        ${exclude_bed} \
        ${exclude_padding} \
        ${output_invariant}

    """
}