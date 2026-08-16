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

    Rscript "${projectDir}/bin/plot_pca.R" \
        "${projectDir}" \
        "${params.rdata}" \
        "${eigval}" \
        "${eigvec}" \
        "${popmap}" 

    """
}