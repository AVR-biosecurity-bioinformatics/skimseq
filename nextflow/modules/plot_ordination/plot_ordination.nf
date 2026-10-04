process PLOT_ORDINATION {
    tag "${distmat}"
    conda "${moduleDir}/environment.yml"

    input:
    path(distmat)
    path(popmap)
    val(covariance)

    output: 
    path("*.pdf"),             emit: plots

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
            library(ggplot2)
        })


        prefix <- "${distmat.simpleName}"
        plot_file <- paste0(prefix, "_mds.pdf")

        M <- read_distance_matrix("${distmat}")

        popmap <- read_popmap("${popmap}")

        if (nrow(M) >= 3L) {
            mds <- cmdscale(
                as.dist(M),
                k = min(2L, nrow(M) - 1L),
                eig = TRUE,
                add = FALSE
            )

            coordinates <- as.data.frame(mds\$points)

            # Degenerate matrices may return only one coordinate.
            if (ncol(coordinates) == 1L) {
                coordinates[[2]] <- 0
            }

            coordinates <- coordinates[, 1:2, drop = FALSE]
            names(coordinates) <- c("PC1", "PC2")

            plot_data <- data.frame(
                sample = rownames(coordinates),
                PC1 = coordinates\$PC1,
                PC2 = coordinates\$PC2,
                stringsAsFactors = FALSE
            )

            plot_data\$pop <- popmap\$pop[match(plot_data\$sample, popmap\$sample)]

            positive_eigenvalues <- mds\$eig[mds\$eig > 0]

            variance_explained <- if (length(positive_eigenvalues) > 0L) {
                100 * positive_eigenvalues / sum(positive_eigenvalues)
            } else {
                numeric()
            }

            axis_label <- function(axis, index) {
                percentage <- if (length(variance_explained) >= index) {
                    sprintf("%.1f%%", variance_explained[[index]])
                } else {
                    "NA"
                }
                paste0(axis, " (", percentage, ")")
            }

            ordination_plot <- ggplot(plot_data, aes(x = -PC1, y = PC2, colour = pop )) +
                geom_point(size = 2) +
                labs(
                    x = axis_label("PC1", 1L),
                    y = axis_label("PC2", 2L),
                    colour = "Population"
                ) +
                theme_classic() +
                theme(
                    legend.position = "right"
                )
                
        } else {
            message("Insufficient samples for MDS after filtering: ", nrow(M))

            ordination_plot <- ggplot() +
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
            filename = plot_file,
            plot = ordination_plot,
            width = 11,
            height = 8,
            units = "in"
        )
    RSCRIPT
    """
}