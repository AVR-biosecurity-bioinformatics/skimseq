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
      filter(!is.na(BIN)) %>%
      dplyr::group_by(RULE, FILTER, VARIANT_TYPE, BIN) %>% # Sum the bins from multiple parallel chunks
      dplyr::summarise(COUNT = sum(COUNT)) %>%
      ungroup() %>%
      group_by(VARIANT_TYPE, RULE) %>%
      mutate(PROP = COUNT / sum(COUNT, na.rm = TRUE)) %>%
      mutate(n_grp = dplyr::n()) %>% # Count number of groups per facet
      ungroup()

    # Make sure every facet has more than one record to avoid bug with scale_x_binned
    df_aug <- df %>%
      bind_rows(
        df %>%
          group_by(VARIANT_TYPE, RULE) %>%
          filter(dplyr::n() == 1) %>% # singleton facets only
          slice(1) %>% # duplicate a single row
          ungroup() %>%
          mutate(
            BIN = BIN + 1e-9,
            PROP = 0 # zero height so it won't change the plot
          )
      )

    variant_types <- factor(
      unique(df$VARIANT_TYPE),
      levels = c("snp", "indel", "invariant")
    )
    variant_qc_plots <- vector("list", length = length(variant_types))
    for (v in 1:length(variant_types)) {
      variant_type <- variant_types[v]
      variant_qc_plots[[v]] <- df_aug %>%
        filter(VARIANT_TYPE == variant_type) %>%
        ggplot(aes(x = BIN, y = PROP, fill = FILTER)) +
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
        facet_wrap(VARIANT_TYPE ~ RULE, scales = "free") +
        scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
        scale_y_continuous(labels = scales::percent) +
        scale_x_binned(n.breaks = 50, breaks = scales::breaks_pretty(n = 6)) +
        theme_classic() +
        labs(x = NULL, y = "Proportion") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    }
    # Write out plots
    pdf("variant_filter_qc.pdf", width = 11, height = 8)
    purrr::walk(variant_qc_plots, plot)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_SITE_FILTERS.rda")
    }
  }
)
