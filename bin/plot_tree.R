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

    # Mat is square (samples x samples)
    drop <- rowSums(!is.na(M)) == 1
    M_clean <- M[!drop, !drop, drop = FALSE]

    M_clean <- M[!drop, !drop, drop = FALSE]

    # Check if matrix still has dimensions
    if (!any(dim(M_clean) == 0)) {
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
          label = "All comparisons are NAN",
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
