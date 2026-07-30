process PLINK_REL {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_rel", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname), path("*.rel*"),                           emit: rel
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Then create GRM relationship matrix
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile ${outname} \
        --allow-extra-chr \
        --make-rel square \
        --out ${outname}
    """
}