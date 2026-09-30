process PLOT_PCA {

    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(eigval), path(eigvec)
    path popmap

    output:
    tuple val(outname),
          path("${outname}_pca.pdf"),
          emit: plots

    tuple val(outname),
          path("${outname}_pca.tsv"),
          emit: coordinates

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    functions_r=\$(command -v functions.R) || {
        echo "ERROR: functions.R not found on PATH" >&2
        exit 1
    }

    Rscript --vanilla - "\$functions_r" <<'RSCRIPT'
        suppressPackageStartupMessages({
            library(ggplot2)
            library(readr)
        })

        functions_file <- commandArgs(trailingOnly = TRUE)[[1]]
        source(functions_file, local = TRUE)

        # Read and validate the PLINK PCA outputs.
        pca <- read_plink_pca(
            eigenvec_file = "${eigvec}",
            eigenval_file = "${eigval}"
        )

        plot_data <- pca\$coordinates
        eigenvalues <- pca\$eigenvalues

        # Calculate variance explained using positive eigenvalues.
        positive_eigenvalues <- pmax(eigenvalues, 0)
        total_variance <- sum(positive_eigenvalues)

        variance_explained <- if (is.finite(total_variance) && total_variance > 0 ) {
            100 * positive_eigenvalues / total_variance
        } else {
            rep(NA_real_, length(eigenvalues))
        }

        pc_percent <- sprintf(
            "%.1f%%",
            variance_explained[1:2]
        )

        pc_percent[ !is.finite(variance_explained[1:2]) ] <- "NA"

        # Read and validate the population map.
        popmap <- read_popmap("${popmap}")

        # Add population assignments to the PCA coordinates.
        plot_data\$pop <- popmap\$pop[match(plot_data\$sample, popmap\$sample)]

        plot_data\$pop[is.na(plot_data\$pop) | plot_data\$pop == "" ] <- "Unknown"

        if (nrow(plot_data) >= 2L) {
            pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, colour = pop )) +
                geom_point(size = 2) +
                labs(
                    x = paste0(
                        "PC1 (",
                        pc_percent[[1]],
                        ")"
                    ),
                    y = paste0(
                        "PC2 (",
                        pc_percent[[2]],
                        ")"
                    ),
                    colour = "Population"
                ) +
                theme_classic() +
                theme(
                    legend.position = "right"
                )
        } else {
            message(
                "Insufficient samples for PCA plot: ",
                nrow(plot_data)
            )

            pca_plot <- ggplot() +
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
            filename = "${outname}_pca.pdf",
            plot = pca_plot,
            device = "pdf",
            width = 11,
            height = 8,
            units = "in"
        )

        readr::write_tsv(
            plot_data,
            "${outname}_pca.tsv"
        )
    RSCRIPT
    """
}