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
    df <- readr::read_tsv(
      table_files,
      col_names = c("RULE", "FILTER", "VARIANT_TYPE", "BIN", "COUNT"),
      col_types = c("cccnn")
    ) %>%
      dplyr::group_by(RULE, FILTER, VARIANT_TYPE, BIN) %>%
      dplyr::summarise(COUNT = sum(COUNT))

    gg.filter_qc <- df %>%
      ggplot(aes(x = BIN, y = COUNT, fill = FILTER)) +
      geom_col() +
      #TODO: add back in these vlines for filter thresholds
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
      facet_wrap(VARIANT_TYPE ~ RULE, scales = "free_x") +
      scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none")

    # Write out plots
    pdf("variant_filter_qc.pdf", width = 11, height = 8)
    plot(gg.filter_qc)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_SITE_FILTERS.rda")
    }
  }
)
