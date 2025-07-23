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
      "vegan",
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
    popmap <- read.table(
      popmap_file,
      header = FALSE,
      col.names = c("sample", "pop")
    )

    # Mat is square (samples x samples)
    drop <- apply(
      M,
      1,
      function(r, i) {
        # remove the self-comparison (diagonal) before testing
        all(is.na(r[-i]) | is.nan(r[-i]))
      },
      i = seq_len(nrow(M))
    )

    M_clean <- M[!drop, !drop, drop = FALSE]

    # Check if matrix still has dimensions
    if (!any(dim(M_clean) == 0)) {
      # Convert matrix to dist matrix
      distmat <- as.dist(M_clean)

      # Conduct MDS
      mds <- capscale(distmat ~ 1)

      # Get target principal components to plot
      target_pc <- c("PC1", "PC2")

      # Create variance explained labels
      var_exp <- scales::percent(mds$CA$eig / sum(mds$CA$eig), accuracy = 0.1)
      names(var_exp) <- names(var_exp) %>% str_replace("MDS", "PC")
      lab_x <- paste0(
        target_pc[1],
        " (",
        var_exp[match(target_pc[1], names(var_exp))],
        ")"
      )
      lab_y <- paste0(
        target_pc[2],
        " (",
        var_exp[match(target_pc[2], names(var_exp))],
        ")"
      )

      # Extract data for ordination
      pcx <- as.data.frame(mds$Ybar)
      colnames(pcx) <- paste0("PC", seq(1, ncol(pcx), 1))

      # Plot ordination
      gg.ord <- pcx %>%
        dplyr::select(paste0("PC", seq(1, 2, 1))) %>% #select top 2 PCs
        tibble::rownames_to_column("sample") %>%
        dplyr::left_join(popmap) %>%
        ggplot(aes(-get(target_pc[1]), get(target_pc[2]), colour = pop)) +
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
          label = "All comparisons are NAN",
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
