#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]
    mat_file <- args[3]
    popmap_file <- args[4]
    params.covariance <- args[5]

    sys.source(paste0(projectDir, "/bin/functions.R"), envir = .GlobalEnv)

    ### load only required packages
    process_packages <- c(
      "dplyr",
      "tidyr",
      "readr",
      "ggplot2",
      "stringr",
      NULL
    )
    invisible(lapply(
      process_packages,
      library,
      character.only = TRUE,
      warn.conflicts = FALSE
    ))

    ### run code

    # Get prefix for output
    prefix <- mat_file %>% stringr::str_remove("\\..*$")

    # Read in matrix files
    M <- read.table(mat_file, row.names = 1)
    colnames(M) <- rownames(M)

    # Read in popmap file
    popmap <- readr::read_tsv(
      popmap_file,
      col_names = c("sample", "pop"),
      col_types = c("cc")
    )

    # Create function that drops sample with most missing data until matrix is complete
    drop_until_complete <- function(M, verbose = TRUE) {
      M <- as.matrix(M)
      M[is.nan(M)] <- NA # normalise NaN to NA

      dropped <- character(0)

      repeat {
        # how many NA per row/col
        na_r <- rowSums(is.na(M))
        na_c <- colSums(is.na(M))

        # max NA across rows/cols
        max_na <- max(c(na_r, na_c), na.rm = TRUE)

        # Break if no NA left
        if (max_na == 0) {
          break
        }

        # pick the worst sample
        worst <- names(which.max(na_r))

        if (verbose) {
          message("Dropping: ", worst, " (", max_na, " NA)")
        }

        # drop that row & col
        M <- M[
          setdiff(rownames(M), worst),
          setdiff(colnames(M), worst),
          drop = FALSE
        ]
        dropped <- c(dropped, worst)

        # safety: if we fall below 2 samples, stop
        if (nrow(M) < 2) break
      }

      # clean up: enforce symmetry & 0 diagonal
      if (nrow(M) > 1) {
        M[is.na(M)] <- 0
        M <- (M + t(M)) / 2
        diag(M) <- 0
      }

      return(M)
    }

    M_clean <- drop_until_complete(M)

    # Check if matrix has enough dimensions
    if (all(dim(M_clean) > 2)) {
      # Convert matrix to dist matrix
      distmat <- as.dist(M_clean)

      # Classical metric MDS / principal coordinates analysis.
      mds <- stats::cmdscale(
        distmat,
        k = min(2L, nrow(M_clean) - 1L),
        eig = TRUE,
        add = FALSE
      )

      pcx <- as.data.frame(mds$points)

      # cmdscale can return only one usable dimension in degenerate cases.
      if (ncol(pcx) == 1) {
        pcx$Dim2 <- 0
      }

      pcx <- pcx[, seq_len(min(2L, ncol(pcx))), drop = FALSE]
      colnames(pcx) <- c("PC1", "PC2")[seq_len(ncol(pcx))]

      if (!"PC2" %in% colnames(pcx)) {
        pcx$PC2 <- 0
      }

      pcx <- tibble::rownames_to_column(pcx, "sample")

      # Calculate percentages from positive eigenvalues only.
      positive_eigenvalues <- mds$eig[mds$eig > 0]

      if (length(positive_eigenvalues) > 0) {
        variance_explained <- 100 * positive_eigenvalues /
          sum(positive_eigenvalues)
      } else {
        variance_explained <- numeric()
      }

      pc1_percent <- if (length(variance_explained) >= 1) {
        sprintf("%.1f%%", variance_explained[1])
      } else {
        "NA"
      }

      pc2_percent <- if (length(variance_explained) >= 2) {
        sprintf("%.1f%%", variance_explained[2])
      } else {
        "NA"
      }

      lab_x <- paste0("PC1 (", pc1_percent, ")")
      lab_y <- paste0("PC2 (", pc2_percent, ")")

      plot_data <- pcx %>%
        dplyr::left_join(popmap, by = "sample")

      # Plot ordination
      gg.ord <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = -.data$PC1, y = .data$PC2, colour = .data$pop)
        ) +
        geom_point(size = 2) +
        labs(x = lab_x, y = lab_y) +
        theme_classic() +
        theme(legend.position = "right")
    } else {
      # If all are dropped by NAN filter, create empty plot
      gg.ord <- ggplot() +
        xlim(0, 1) +
        ylim(0, 1) + # give coords to place the text
        annotate(
          "text",
          x = .5,
          y = .5,
          label = "Insufficient samples to make plot",
          size = 6,
          fontface = "bold"
        ) +
        theme_void()
    }

    # Write out plots
    pdf(paste0(prefix, "_ord.pdf"), width = 11, height = 8)
    plot(gg.ord)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_ORDINATION.rda")
    }
  }
)
