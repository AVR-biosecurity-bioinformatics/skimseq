process JOINT_GENOTYPE {
    def process_name = "joint_genotype"    
    // tag "-"
    // label "small"
    time '2.h'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_list), path(gvcf), path(gvcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${gvcf} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_list}

    """
}