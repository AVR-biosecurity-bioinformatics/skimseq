process PLINK_REL {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_rel", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname),
          path("${outname}.rel"),
          path("${outname}.rel.id"),
          emit: rel

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMPLES=\$(wc -l < "${outname}.fam")
    N_VARIANTS=\$(wc -l < "${outname}.bim")

    if (( N_SAMPLES < 2 )); then
        echo "ERROR: Relationship matrix requires at least two samples" >&2
        exit 1
    fi

    if (( N_VARIANTS < 1 )); then
        echo "ERROR: Relationship matrix requires at least one variant" >&2
        exit 1
    fi

    # Generate explicit allele frequencies. This permits small test datasets
    # to run and records the frequencies used to standardise the GRM.
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --freq \
        --out "${outname}"

    # Generate the variance-standardised genomic relationship matrix.
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --read-freq "${outname}.afreq" \
        --make-rel square \
        --out "${outname}"
    """
}