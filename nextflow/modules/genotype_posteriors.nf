process GENOTYPE_POSTERIORS {
    publishDir "${launchDir}/output/modules/genotype_posteriors", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(tbi), path(sites_vcf), path(sites_tbi)


    output: 
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path("*.gp.vcf.gz"), path("*.gp.vcf.gz.tbi"),      emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash genotype_posteriors.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} \
        ${vcf} \
        ${sites_vcf}

    """
}