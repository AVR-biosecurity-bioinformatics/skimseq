process SPLIT_VCF_BY_TYPE {
    tag "${outname}"
    publishDir "${launchDir}/output/modules/split_vcf_by_type", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),
          path("${outname}.snp.vcf.gz"),
          path("${outname}.snp.vcf.gz.tbi"),
          emit: snp_vcf

    tuple val(outname),
          path("${outname}.indel.vcf.gz"),
          path("${outname}.indel.vcf.gz.tbi"),
          emit: indel_vcf

    tuple val(outname),
          path("${outname}.invariant.vcf.gz"),
          path("${outname}.invariant.vcf.gz.tbi"),
          emit: invariant_vcf

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    bcftools view -Ou ${vcf} \
    | tee \
        >(bcftools view -Oz -v snps   -o ${outname}.snp.vcf.gz) \
        >(bcftools view -Oz -v indels -o ${outname}.indel.vcf.gz) \
    | bcftools view -Oz -v ref -o ${outname}.invariant.vcf.gz

    # Index outputs
    bcftools index -t --threads ${task.cpus} ${outname}.snp.vcf.gz
    bcftools index -t --threads ${task.cpus} ${outname}.indel.vcf.gz
    bcftools index -t --threads ${task.cpus} ${outname}.invariant.vcf.gz
    """
}