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

    # List matrix files
    mat_files <- list.files(pattern = "\\.mat$")

    prefix <- mat_files %>% stringr::str_remove("\\..*$")

    # Read in matrix files
    mat <- read.table(mat_files, row.names = 1)
    colnames(mat) <- rownames(mat)

    # Convert to distance matrix
    distmat <- as.dist(mat)

    # TODO: Check if the input is a convariance matrix, if so run eigen function

    # TODO: Handle population labels - if missing colour by sample

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
      dplyr::select(paste0("PC", seq(1, 2, 1))) %>%
      tibble::rownames_to_column("sample") %>%
      ggplot(aes(-get(target_pc[1]), get(target_pc[2]))) + #, colour = pop
      geom_point(size = 2) +
      labs(x = lab_x, y = lab_y) +
      theme_classic() +
      theme(legend.position = "none")

    # Write out plots
    pdf(paste0(prefix, "_ord.pdf"), width = 11, height = 8)
    plot(gg.ord)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "CREATE_INTERVALS.rda")
    }
  }
)
