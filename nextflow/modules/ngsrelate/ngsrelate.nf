process NGSRELATE {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/ngsrelate", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/ngsrelate", mode: 'copy'

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
        -O ${prefix}.res \
        -I 1 
    """
}