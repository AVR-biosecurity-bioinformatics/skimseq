process GENOTYPE_POSTERIORS {
    def process_name = "genotype_posteriors"    
    // tag "-"
    // label "small"
    time '2.h'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_list), path(gvcf), path(gvcf_tbi)


    output: 
    tuple val(interval_hash), path(interval_list), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),      emit: gvcf_intervals
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${gvcf} \
        ${interval_hash} \
        ${interval_list}

    """
}