nextflow.enable.types = true

process PLINK_PCA {

    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        bed: Path,
        bim: Path,
        fam: Path
    )

    output:
    record(
        outname: outname,
        eigenval: file("${outname}.eigenval"),
        eigenvec: file("${outname}.eigenvec")
    )

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
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --freq \
        --out "${outname}"

    # Run PLINK pca
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --read-freq "${outname}.afreq" \
        --pca "\${N_PCS}" \
        --out "${outname}"

    """
}