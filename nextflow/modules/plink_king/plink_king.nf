nextflow.enable.types = true

process PLINK_KING {

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
        king: files("${outname}.king*")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create KING relationship table.
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --bfile "${outname}" \
        --allow-extra-chr \
        --make-king-table \
        --out "${outname}.king"
    """
}