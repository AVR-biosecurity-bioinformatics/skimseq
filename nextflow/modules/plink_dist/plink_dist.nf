nextflow.enable.types = true

process PLINK_DIST {

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
        mat: file("${outname}.mat")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMPLES=\$(wc -l < "${fam}")
    N_VARIANTS=\$(wc -l < "${bim}")

    if (( N_SAMPLES < 2 )); then
        echo "ERROR: distance matrix requires at least two samples" >&2
        exit 1
    fi

    if (( N_VARIANTS < 1 )); then
        echo "ERROR: distance matrix requires at least one variant" >&2
        exit 1
    fi

    plink \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --distance square 1-ibs flat-missing \
        --out tmp

    paste \
        <(awk '{ print \$2 }' tmp.mdist.id) \
        tmp.mdist \
        > "${outname}.mat"

    rm -f \
        tmp.mdist \
        tmp.mdist.id \
        tmp.log \
        tmp.nosex
    """
}