process PLOT_TREE {
    def process_name = "plot_tree"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(distmat)

    output: 
    path("*.pdf"),             emit: plots

    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata}
    """
}