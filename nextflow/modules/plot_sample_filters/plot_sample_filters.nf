process PLOT_SAMPLE_FILTERS {
    publishDir "${launchDir}/output/modules/plot_sample_filters", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'

    input:
    path(missing_summary)
    val(sample_max_missing)

    output: 
    path("*.pdf"),               emit: plots
    path("sample_missing.tsv"),  emit: tsv

    script:
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        Rscript ${projectDir}/bin/plot_sample_filters.R \
        ${projectDir} \
        ${params.rdata} \
        "${sample_max_missing}" 
    """
}