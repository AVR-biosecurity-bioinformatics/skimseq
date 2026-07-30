process PLINK_IMPORT {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_import", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    tuple val(outname), path("${outname}.{bim,bed,fam}"),                           emit: plink
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create PLINK bed file
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --vcf ${vcf} \
        --allow-extra-chr \
        --double-id \
        --make-bed \
        --out ${outname}

    """
}