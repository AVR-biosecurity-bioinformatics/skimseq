process PLOT_ORDINATION {
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_ordination", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(distmat)
    path(popmap)
    val(covariance)

    output: 
    path("*.pdf"),             emit: plots

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    Rscript "${projectDir}/bin/plot_ordination.R" \
        "${projectDir}" \
        "${params.rdata}" \
        "${distmat}" \
        "${popmap}" \
        "${covariance}"
    """
}