#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]
    mat_file <- args[3]
    popmap_file <- args[4]

    # TODO: better to make some kind of temporary popmap file, potentially using nextflow native commands

    sys.source(paste0(projectDir, "/bin/functions.R"), envir = .GlobalEnv)

    ### load only required packages
    process_packages <- c(
      "tibble",
      "dplyr",
      "ape",
      "ggtree",
      "ggplot2",
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

      # Construct a tree using NJ method
      tree <- ape::nj(distmat)

      # save the tree to Newick format file
      write.tree(tree, file = paste0(prefix, "_tree.nwk"))

      # Create base tree
      p1 <- ggtree(tree, layout = "equal_angle") # see more at ggtree

      # Add population colours
      tree_df <- tibble::enframe(
        tree$tip.label,
        name = NULL,
        value = "sample"
      ) %>%
        dplyr::left_join(popmap)

      # Update tree with colours
      gg.tree <- p1 %<+% (tree_df) + aes(color = pop)
    } else {
      # If all are dropped by NAN filter, create empty plot
      gg.tree <- ggplot() +
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
    pdf(paste0(prefix, "_tree.pdf"), width = 11, height = 8)
    plot(gg.tree)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_TREE.rda")
    }
  }
)
