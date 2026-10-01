process NGSRELATE {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(bcf), path(csi)

    output: 
    tuple val(outname),
          path("${outname}.rel"),
          path("${outname}.rel.id"),
          emit: rel

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