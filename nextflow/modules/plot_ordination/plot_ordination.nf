process PLOT_ORDINATION {
    tag "${distmat}"
    conda "${moduleDir}/environment.yml"

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

    plot_ordination.R \
        "${params.rdata}" \
        "${distmat}" \
        "${popmap}" \
        "${covariance}"
    """
}