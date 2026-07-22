process SPLIT_VCF_BY_TYPE {
    def process_name = "split_vcf_by_type"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.snp.{vcf,g.vcf}.gz"), path("${outname}.snp.{vcf,g.vcf}.gz.tbi"),       emit: snp_vcf
    tuple val(outname),  path("${outname}.indel.{vcf,g.vcf}.gz"), path("${outname}.indel.{vcf,g.vcf}.gz.tbi"),       emit: indel_vcf
    tuple val(outname),  path("${outname}.invariant.{vcf,g.vcf}.gz"), path("${outname}.indel.{vcf,g.vcf}.gz.tbi"),       emit: invariant_vcf

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${outname}"

    """
}