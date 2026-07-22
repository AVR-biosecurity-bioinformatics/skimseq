process SPLIT_VCF_BY_TYPE {
    publishDir "${launchDir}/output/modules/split_vcf_by_type", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.snp.{vcf,g.vcf}.gz"), path("${outname}.snp.{vcf,g.vcf}.gz.tbi"),       emit: snp_vcf
    tuple val(outname),  path("${outname}.indel.{vcf,g.vcf}.gz"), path("${outname}.indel.{vcf,g.vcf}.gz.tbi"),       emit: indel_vcf
    tuple val(outname),  path("${outname}.invariant.{vcf,g.vcf}.gz"), path("${outname}.indel.{vcf,g.vcf}.gz.tbi"),       emit: invariant_vcf

    script:
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash split_vcf_by_type.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${outname}"

    """
}