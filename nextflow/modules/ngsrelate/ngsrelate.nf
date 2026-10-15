nextflow.enable.types = true

process NGSRELATE {

    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        bcf: Path,
        csi: Path
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

    ngsRelate \
        -p ${task.cpus} \
        -h ${bcf} \
        -O ${outname}.res \
        -I 1 
    """
}