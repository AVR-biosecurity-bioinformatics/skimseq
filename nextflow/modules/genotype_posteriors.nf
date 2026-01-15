process GENOTYPE_POSTERIORS {
    def process_name = "genotype_posteriors"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:BCFtools/1.22-GCC-13.3.0"

    input:
    tuple val(interval_hash), path(interval_bed), path(vcf), path(tbi), path(sites_vcf), path(sites_tbi)


    output: 
    tuple val(interval_hash), path(interval_bed), path("*.vcf.gz"), path("*.vcf.gz.tbi"),      emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} \
        ${vcf} \
        ${sites_vcf}

    """
}