process SPLIT_VCF_BY_TYPE {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(bcf), path(csi)
    
    output: 
    tuple val(outname),
          path("${outname}.snp.bcf"),
          path("${outname}.snp.bcf.csi"),
          emit: snp_vcf

    tuple val(outname),
          path("${outname}.indel.bcf"),
          path("${outname}.indel.bcf.csi"),
          emit: indel_vcf

    tuple val(outname),
          path("${outname}.invariant.bcf"),
          path("${outname}.invariant.bcf.csi"),
          emit: invariant_vcf

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    bcftools view -Ou ${bcf} \
    | tee \
        >(bcftools view -Oz -v snps   -o ${outname}.snp.bcf) \
        >(bcftools view -Oz -v indels -o ${outname}.indel.bcf) \
    | bcftools view -Ob -v ref -o ${outname}.invariant.bcf

    # Index outputs
    bcftools index --threads ${task.cpus} ${outname}.snp.bcf
    bcftools index --threads ${task.cpus} ${outname}.indel.bcf
    bcftools index --threads ${task.cpus} ${outname}.invariant.bcf
    """
}