#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]

    sys.source(paste0(projectDir, "/bin/functions.R"), envir = .GlobalEnv)

    ### load only required packages
    process_packages <- c(
      "dplyr",
      "tidyr",
      "readr",
      "ggplot2",
      "stringr",
      "patchwork",
      NULL
    )
    invisible(lapply(
      process_packages,
      library,
      character.only = TRUE,
      warn.conflicts = FALSE
    ))

    ### run code

    # List filtering summary files
    table_files <- list.files(pattern = "\\.tsv.gz$")

    # Read in all tables and combine
    df <- read_tsv(
      table_files,
      col_names = c("name", "label", "type", "bin", "count")
    )

    gg.site_qc_dist <- df %>%
      ggplot(aes(x = bin, y = count, fill = label)) +
      geom_col() +
      #geom_density()+
      #geom_vline(
      #  data = parameter_table %>% dplyr::slice(i),
      #  aes(xintercept = lower),
      #  lty = "dashed"
      #) +
      #geom_vline(
      #  data = parameter_table %>% dplyr::slice(i),
      #  aes(xintercept = upper),
      #  lty = "dashed"
      #) +
      facet_wrap(type ~ name) +
      scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
      labs(
        x = paste0(type_to_select, " ", filter_to_select),
        y = "Count"
      ) +
      theme_classic() +
      theme(legend.position = "none")

    # Write out plots
    pdf("site_filtering_qc.pdf", width = 11, height = 8)
    plot(gg.site_qc_dist)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_SITE_FILTERS.rda")
    }
  }
)
