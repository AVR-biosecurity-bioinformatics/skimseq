nextflow.enable.types = true

process PLINK_IMPORT {

    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        vcf: Path,
        index: Path
    )

    output:
    record(
        outname: outname,
        bed: file("${outname}.bed"),
        bim: file("${outname}.bim"),
        fam: file("${outname}.fam")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create PLINK bed file
    # TODO: Add filters
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --bcf ${vcf} \
        --mind 0.9 \
        --allow-extra-chr \
        --double-id \
        --make-bed \
        --out ${outname}

    """
}