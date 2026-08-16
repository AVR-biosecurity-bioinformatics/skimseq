process PLOT_TREE {
    tag "${distmat}"
    conda "${moduleDir}/environment.yml"

    input:
    path(distmat)
    path(popmap)

    output: 
    path("*.pdf"),             emit: plots
    path("*.nwk"),             emit: newick_tree

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    plot_tree.R \
        ${params.rdata} \
        ${distmat} \
        ${popmap}
    """
}