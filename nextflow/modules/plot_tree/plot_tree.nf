process PLOT_TREE {
    //conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/plot_tree", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/visualisation/trees", mode: 'copy'

    input:
    path(distmat)
    path(popmap)

    output: 
    path("*.pdf"),             emit: plots

    script:
    """
    shifter --image=gmboowa/ggtree:latest -- \
        Rscript ${projectDir}/bin/plot_tree.R \
        ${projectDir} \
        ${params.rdata} \
        ${distmat} \
        ${popmap}
    """
}