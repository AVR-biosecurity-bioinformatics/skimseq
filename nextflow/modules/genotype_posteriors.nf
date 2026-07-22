process GENOTYPE_POSTERIORS {
    def process_name = "genotype_posteriors"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(tbi), path(sites_vcf), path(sites_tbi)


    output: 
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path("*.gp.vcf.gz"), path("*.gp.vcf.gz.tbi"),      emit: vcf
    
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