process PLOT_SAMPLE_FILTERS {
    tag "${missing_summary}"
    conda "${moduleDir}/environment.yml"

    input:
    path(missing_summary)

    output: 
    path("*.pdf"),               emit: plots
    path("sample_missing.tsv"),  emit: sample_missing_tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    Rscript ${projectDir}/bin/plot_sample_filters.R \
        ${params.rdata} \
        "${params.vcf_sample_max_missing}" 
    """
}