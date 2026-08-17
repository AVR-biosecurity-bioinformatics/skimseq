process NGSRELATE {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

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
        -h ${vcf} \
        -O ${outname}.res \
        -I 1 
    """
}