#!/usr/bin/env Rscript
.libPaths(.Library)

params.rdata <- "false"

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    params.rdata <- args[1]
    eigenval_file <- args[2]
    eigenvec_file <- args[3]
    popmap_file <- args[4]


    process_packages <- c(
      "dplyr",
      "readr",
      "ggplot2",
      "stringr"
    )

    invisible(lapply(
      process_packages,
      library,
      character.only = TRUE,
      warn.conflicts = FALSE
    ))

    prefix <- eigenvec_file %>%
      basename() %>%
      stringr::str_remove("\\.eigenvec(?:\\.zst)?$")

    eigenvec <- readr::read_table(
      eigenvec_file,
      show_col_types = FALSE,
      progress = FALSE
    )

    # PLINK 2 normally writes #FID, IID, and PC columns. Remove the leading
    # hash from the first field name and use IID as the sample identifier.
    names(eigenvec) <- stringr::str_remove(names(eigenvec), "^#")

    if (!"IID" %in% names(eigenvec)) {
      stop(
        "The PLINK eigenvector file does not contain an IID column: ",
        eigenvec_file
      )
    }

    pc_columns <- grep(
      "^PC[0-9]+$",
      names(eigenvec),
      value = TRUE
    )

    if (length(pc_columns) < 2) {
      stop(
        "At least two principal components are required, but found: ",
        paste(pc_columns, collapse = ", ")
      )
    }

    eigenvalues <- readr::read_lines(
      eigenval_file,
      progress = FALSE
    ) %>%
      trimws() %>%
      as.numeric()

    if (
      length(eigenvalues) < 2 ||
      anyNA(eigenvalues[1:2])
    ) {
      stop(
        "The PLINK eigenvalue file does not contain two valid eigenvalues: ",
        eigenval_file
      )
    }

    positive_eigenvalues <- pmax(eigenvalues, 0)
    total_variance <- sum(positive_eigenvalues)

    if (is.finite(total_variance) && total_variance > 0) {
      variance_explained <- 100 * positive_eigenvalues / total_variance
    } else {
      variance_explained <- rep(NA_real_, length(eigenvalues))
    }

    pc1_percent <- if (is.finite(variance_explained[1])) {
      sprintf("%.1f%%", variance_explained[1])
    } else {
      "NA"
    }

    pc2_percent <- if (is.finite(variance_explained[2])) {
      sprintf("%.1f%%", variance_explained[2])
    } else {
      "NA"
    }

    popmap <- readr::read_tsv(
      popmap_file,
      col_names = c("sample", "pop"),
      col_types = "cc",
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      dplyr::distinct(.data$sample, .keep_all = TRUE)

    duplicated_samples <- eigenvec$IID[
      duplicated(eigenvec$IID)
    ]

    if (length(duplicated_samples) > 0) {
      stop(
        "Duplicate sample IDs in PLINK eigenvector file: ",
        paste(unique(duplicated_samples), collapse = ", ")
      )
    }

    plot_data <- eigenvec %>%
      dplyr::transmute(
        sample = as.character(.data$IID),
        PC1 = as.numeric(.data$PC1),
        PC2 = as.numeric(.data$PC2)
      ) %>%
      dplyr::left_join(popmap, by = "sample") %>%
      dplyr::mutate(
        pop = dplyr::coalesce(.data$pop, "Unknown")
      )

    if (nrow(plot_data) >= 2) {
      gg.pca <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$PC1,
          y = .data$PC2,
          colour = .data$pop
        )
      ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(
          x = paste0("PC1 (", pc1_percent, ")"),
          y = paste0("PC2 (", pc2_percent, ")"),
          colour = "Population"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = "right"
        )
    } else {
      gg.pca <- ggplot2::ggplot() +
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "Insufficient samples to make plot",
          size = 6,
          fontface = "bold"
        ) +
        ggplot2::theme_void()
    }

    ggplot2::ggsave(
      filename = paste0(prefix, "_pca.pdf"),
      plot = gg.pca,
      width = 11,
      height = 8,
      units = "in"
    )

    readr::write_tsv(
      plot_data,
      paste0(prefix, "_pca.tsv")
    )
  },
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1)
  },
  finally = {
    if (identical(params.rdata, "true")) {
      save.image(file = "PLOT_PLINK_PCA.rda")
    }
  }
)