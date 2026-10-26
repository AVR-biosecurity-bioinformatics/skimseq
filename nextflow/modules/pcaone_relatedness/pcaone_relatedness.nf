nextflow.enable.types = true

process PCAONE_RELATEDNESS {

    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        bed: Path,
        bim: Path,
        fam: Path
    )
    relatedness_n_pcs: Integer

    output:
    record(
        outname: outname,
        relatedness: file("${outname}.relatedness.kinship"),
        correlation: file("${outname}.relatedness.corres")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMPLES=\$(wc -l < "${fam}")
    N_VARIANTS=\$(wc -l < "${bim}")

    if (( N_SAMPLES < 2 )); then
        echo "ERROR: Relatedness analysis requires at least two samples" >&2
        exit 1
    fi

    if (( N_VARIANTS < 1 )); then
        echo "ERROR: Relatedness analysis requires at least one variant" >&2
        exit 1
    fi

    N_PCS=20
    (( N_PCS > N_SAMPLES - 1 )) && N_PCS=\$((N_SAMPLES - 1))
    (( N_PCS > N_VARIANTS )) && N_PCS=\${N_VARIANTS}

    N_RELATEDNESS_PCS=${relatedness_n_pcs}

    if (( N_RELATEDNESS_PCS > N_PCS )); then
        echo "ERROR: Requested \${N_RELATEDNESS_PCS} ancestry PCs, but only \${N_PCS} can be calculated" >&2
        exit 1
    fi

    # Calculate the PCA model used for ancestry adjustment.
    PCAone \
        --bfile "${outname}" \
        -k "\${N_PCS}" \
        --threads ${task.cpus} \
        -o "${outname}.relatedness_pca"

    # Estimate relatedness using K-1 ancestry PCs.
    PCAone \
        --bfile "${outname}" \
        -P "${outname}.relatedness_pca" \
        --evaladmix \
        -k "\${N_RELATEDNESS_PCS}" \
        --threads ${task.cpus} \
        -o "${outname}.relatedness"

    #--evaladmix-ibd \
    """
}