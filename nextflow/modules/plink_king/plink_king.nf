process PLINK_KING {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

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