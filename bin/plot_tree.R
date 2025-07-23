#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]
    params.covariance <- args[3]
    #popmap <- args[3]

    sys.source(paste0(projectDir, "/bin/functions.R"), envir = .GlobalEnv)

    ### load only required packages
    process_packages <- c(
      "dplyr",
      "tidyr",
      "readr",
      "ggplot2",
      "stringr",
      "ape",
      "ggtree",
      NULL
    )
    invisible(lapply(
      process_packages,
      library,
      character.only = TRUE,
      warn.conflicts = FALSE
    ))

    ### run code

    # List matrix files
    mat_files <- list.files(pattern = "\\.mat$")

    prefix <- mat_files %>% stringr::str_remove("\\..*$")

    # Read in matrix files
    mat <- read.table(mat_files, row.names = 1)
    colnames(mat) <- rownames(mat)

    # Convert to dist matrix
    distmat <- as.dist(mat)

    # Construct a tree using NJ method
    tree <- nj(distmat)

    # save the tree to Newick format file
    write.tree(tree, file = paste0(prefix, "_tree.nwk"))

    gg.tree <- ggtree(tree, layout = "equal_angle") # see more at ggtree

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
