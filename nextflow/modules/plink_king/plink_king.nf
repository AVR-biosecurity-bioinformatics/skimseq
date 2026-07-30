process PLINK_KING {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_king", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname), path("*.king*"),                           emit: king
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Then create king table
    plink2 \
    --threads ${task.cpus} \
    --memory ${task.memory.mega} \
    --bfile ${outname} \
    --allow-extra-chr \
    --make-king-table \
    --out ${outname}.king

    """
}