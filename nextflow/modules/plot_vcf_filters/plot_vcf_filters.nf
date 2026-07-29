process PLOT_VCF_FILTERS {
    // tag "${outname}"
    //conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_vcf_filters", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'

    input:
    path(filter_hist)
    //path(filter_summary)
    val(outname)

    output: 
    path("*.pdf"),   emit: plots
    //path("*.tsv"),   emit: summary

    script:
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        Rscript ${projectDir}/bin/plot_vcf_filters.R \
        ${projectDir} \
        ${params.rdata} \
        ${outname}
    """
}