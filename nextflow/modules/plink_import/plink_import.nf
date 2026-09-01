process PLINK_IMPORT {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    tuple val(outname), path("${outname}.{bim,bed,fam}"),                           emit: plink
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create PLINK bed file
    # TODO: Add filters
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --vcf ${vcf} \
        --mind 0.9 \
        --allow-extra-chr \
        --double-id \
        --make-bed \
        --out ${outname}

    """
}