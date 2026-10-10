nextflow.enable.types = true

process PLINK_REL {

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
        rel: file("${outname}.rel"),
        rel_id: file("${outname}.rel.id")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMPLES=\$(wc -l < "${fam}")
    N_VARIANTS=\$(wc -l < "${bim}")

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
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --freq \
        --out "${outname}"

    # Generate the variance-standardised genomic relationship matrix.
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --read-freq "${outname}.afreq" \
        --make-rel square \
        --out "${outname}"
    """
}