process PLOT_ORDINATION {
    def process_name = "plot_ordination"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/visualisation/ordination", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(distmat)
    path(popmap)
    val(covariance)

    output: 
    path("*.pdf"),             emit: plots

    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata} \
        ${distmat} \
        ${popmap} \
        ${covariance}
    """
}