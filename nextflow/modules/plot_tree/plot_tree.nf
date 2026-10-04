process PLOT_TREE {

    tag "${distmat.simpleName}"
    conda "${moduleDir}/environment.yml"

    input:
    path distmat
    path popmap

    output:
    path "${distmat.simpleName}_tree.pdf", emit: plots
    path "${distmat.simpleName}_tree.nwk", emit: newick_tree

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    functions_r=\$(command -v functions.R)

    Rscript --vanilla - "\$functions_r" <<'RSCRIPT'
    
        functions_file <- commandArgs(trailingOnly = TRUE)[[1]]
        source(functions_file, local = TRUE)

        suppressPackageStartupMessages({
            library(ape)
            library(ggplot2)
            library(ggtree)
        })

        prefix <- "${distmat.simpleName}"
        tree_file <- paste0(prefix, "_tree.nwk")
        plot_file <- paste0(prefix, "_tree.pdf")

        M <- read_distance_matrix("${distmat}")

        popmap <- read_popmap("${popmap}")

        if (nrow(M) >= 3L) {
            tree <- ape::nj(as.dist(M))
            ape::write.tree(tree, file = tree_file)

            tree_data <- data.frame(
                sample = tree\$tip.label,
                pop = popmap\$pop[
                    match(tree\$tip.label, popmap\$sample)
                ]
            )

            tree_plot <- ggtree(tree, layout = "equal_angle") %<+%
                tree_data +
                aes(colour = pop)

        } else {
            message(
                "Insufficient samples for NJ tree after filtering: ",
                nrow(M)
            )

            file.create(tree_file)

            tree_plot <- ggplot() +
                annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Insufficient samples to make plot",
                    size = 6,
                    fontface = "bold"
                ) +
                xlim(0, 1) +
                ylim(0, 1) +
                theme_void()
        }

        ggsave(
            plot_file,
            tree_plot,
            width = 11,
            height = 8,
            units = "in"
        )
    RSCRIPT
    """
}