process PLOT_SAMPLE_FILTERS {
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_sample_filters", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(missing_summary)

    output: 
    path("*.pdf"),               emit: plots
    path("sample_missing.tsv"),  emit: tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    Rscript ${projectDir}/bin/plot_sample_filters.R \
        ${projectDir} \
        ${params.rdata} \
        "${params.vcf_sample_max_missing}" 
    """
}