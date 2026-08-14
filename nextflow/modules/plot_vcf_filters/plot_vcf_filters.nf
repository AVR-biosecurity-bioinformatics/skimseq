process PLOT_VCF_FILTERS {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_vcf_filters", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(filter_hist)
    //path(filter_summary)
    val(outname)

    output: 
    path("*.pdf"),   emit: plots
    //path("*.tsv"),   emit: summary

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent loaded HPC Python/R modules from contaminating the Conda environment.
    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    Rscript ${projectDir}/bin/plot_vcf_filters.R \
        ${projectDir} \
        ${params.rdata} \
        ${outname}
    """
}