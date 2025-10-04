process PLOT_SAMPLE_FILTERS {
    def process_name = "plot_sample_filters"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(gt_filter_table)
    val(sample_max_missing)

    output: 
    path("*.pdf"),               emit: plots
    path("sample_missing.tsv"),  emit: tsv

    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata} \
        "${sample_max_missing}" 
    """
}