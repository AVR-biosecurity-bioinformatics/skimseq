#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]

    #snp_filtering parameters
    params.sample_max_missing <- args[3]

    # TESTING
    #params.sample_max_missing <- 0.5

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

    # List table files
    table_files <- list.files(pattern = "\\.tsv$")

    # Read in table file
    df <- read_tsv(table_files) %>%
      dplyr::mutate(
        FILTER = ifelse(
          MISSING_FRACTION >
            as.numeric(params.sample_max_missing %>% str_remove_all("\\[|\\]")),
          "FAIL",
          "PASS"
        )
      )

    gg.sample_qc_dist <- df %>%
      ggplot(aes(x = MISSING_FRACTION, fill = FILTER)) +
      geom_histogram() +
      geom_vline(
        xintercept = as.numeric(
          params.sample_max_missing %>% str_remove_all("\\[|\\]")
        ),
        lty = "dashed"
      ) +
      scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
      labs(
        x = "Missing data",
        y = "Count"
      ) +
      theme_classic() +
      theme(legend.position = "none")

    write_tsv(df, "sample_missing.tsv")

    # Write out plots
    pdf("sample_filtering_qc.pdf", width = 11, height = 8)
    plot(gg.sample_qc_dist)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_SAMPLE_FILTERS.rda")
    }
  }
)
