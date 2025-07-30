process PLOT_GT_FILTERS {
    def process_name = "plot_gt_filters"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(gt_filter_table)
    tuple val(gt_qual), val(gt_dp_min), val(gt_dp_max)

    output: 
    path("*.pdf"),             emit: plots

    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata} \
        "${gt_qual}" \
        "${gt_dp_min}" \
        "${gt_dp_max}" 
    """
}