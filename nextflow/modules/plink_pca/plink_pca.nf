process PLINK_PCA {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plink_pca", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname),
          path("${outname}.pca.eigenval"),
          path("${outname}.pca.eigenvec"),
          emit: pca
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMPLES=\$(wc -l < "${outname}.fam")
    N_VARIANTS=\$(wc -l < "${outname}.bim")

    if (( N_SAMPLES < 2 )); then
        echo "ERROR: PCA requires at least two samples" >&2
        exit 1
    fi

    if (( N_VARIANTS < 1 )); then
        echo "ERROR: PCA requires at least one variant" >&2
        exit 1
    fi

    # PCA rank is limited by the numbers of samples and variants.
    N_PCS=20
    (( N_PCS > N_SAMPLES - 1 )) && N_PCS=\$((N_SAMPLES - 1))
    (( N_PCS > N_VARIANTS )) && N_PCS=\${N_VARIANTS}

    # Generate an explicit frequency file. This permits small test datasets
    # to run and makes the frequencies used by PCA reproducible.
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --freq \
        --out "${outname}.pca"

    # Run PLINK pca
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.mega} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --read-freq "${outname}.pca.afreq" \
        --pca "\${N_PCS}" \
        --out "${outname}.pca"

    """
}