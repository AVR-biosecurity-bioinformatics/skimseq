process PLOT_PCA {
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_pca", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/visualisation/ordination", mode: 'copy'

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