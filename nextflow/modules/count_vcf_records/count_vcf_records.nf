process COUNT_VCF_RECORDS {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/count_vcf_records", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)
    path(interval_bed)
    path(exclude_bed)

    output: 
    tuple val(sample), path("${sample}.counts.bed.gz"),  path("${sample}.counts.bed.gz.tbi"), emit: counts
    tuple val(sample), path("*.missing.tsv"),                                                 emit: missing_frac
    tuple val(sample), path("*dphist.tsv"),                                                   emit: dphist

    script:
    """
    #!/usr/bin/env bash

    ### run process script
    bash count_vcf_records.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${gvcf}" \
        ${ref_genome} \
        ${interval_bed} \
        ${exclude_bed} \
        ${sample} 
        
    """
}
