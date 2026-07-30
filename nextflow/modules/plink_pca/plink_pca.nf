process PLINK_PCA {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_pca", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname), path("*.pca*"),                           emit: pca
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Run PLINK pca
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile ${outname} \
        --allow-extra-chr \
        --pca 20} \
        --out ${chr}.pca

    """
}