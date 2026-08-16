process PLOT_PCA {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(eigval), path(eigvec)
    path(popmap)

    output: 
    path("*.pdf"),             emit: plots

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    plot_pca.R \
        "${params.rdata}" \
        "${eigval}" \
        "${eigvec}" \
        "${popmap}" 

    """
}