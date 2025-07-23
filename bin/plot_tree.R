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
    mat <- read.table(mat_file, row.names = 1)
    colnames(mat) <- rownames(mat)

    # Read in popmap file
    popmap <- read.table(
      popmap_file,
      header = FALSE,
      col.names = c("sample", "pop")
    )

    # Convert matrix to dist matrix
    distmat <- as.dist(mat)

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

    gg.tree <- p1 %<+% (tree_df) + aes(color = pop) + geom_tiplab()

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
