process PLOT_ORDINATION {
    //conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_ordination", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/visualisation/ordination", mode: 'copy'

    input:
    path(distmat)
    path(popmap)
    val(covariance)

    output: 
    path("*.pdf"),             emit: plots

    script:
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        Rscript ${projectDir}/bin/plot_ordination.R \
        ${projectDir} \
        ${params.rdata} \
        ${distmat} \
        ${popmap} \
        ${covariance}
    """
}