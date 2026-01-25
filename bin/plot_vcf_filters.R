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

    hist_files <- list.files(pattern = "hist.tsv.gz$")
    hist_files <- hist_files[file.size(hist_files) > 0]

    # Read in all tables and combine
    df <- readr::read_tsv(
      hist_files,
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

    # Function to not display all breaks for scale_x_binned
    thin_binned_labels <- function(target = 8, fmt = scales::label_number()) {
      force(fmt)
      function(br) {
        n <- length(br)
        if (n <= target) {
          return(fmt(br))
        }
        step <- ceiling(n / target)
        keep <- unique(sort(c(1, seq(1, n, by = step), n))) # first/last + every step
        out <- rep("", n)
        out[keep] <- fmt(br[keep])
        out
      }
    }

    #TODO: add back in  vlines for filter thresholds

    variant_qc_plots <- vector("list", length = length(variant_types))
    for (v in 1:length(variant_types)) {
      variant_type <- variant_types[v]
      variant_qc_plots[[v]] <- df_aug %>%
        filter(VARIANT_TYPE == variant_type) %>%
        ggplot(aes(x = BIN, y = PROP, fill = FILTER)) +
        geom_col() +
        facet_wrap(VARIANT_TYPE ~ RULE, scales = "free") +
        scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
        scale_y_continuous(labels = scales::percent) +
        scale_x_binned(
          n.breaks = 25, # keep many bins
          labels = thin_binned_labels(8) # show ~8 labels per facet
        ) +
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

    # Create a joint table of the summary files
    # List filtering summary files
    summary_files <- list.files(pattern = "filter_summary.tsv$")
    summary_files <- summary_files[file.size(summary_files) > 0]
    df_summary <- readr::read_tsv(
      summary_files,
      col_names = c("FILTER", "COUNT"),
      col_types = c("cn")
    ) %>%
      group_by(FILTER) %>%
      summarise(COUNT = sum(COUNT))

    # Write out summary file
    write_tsv(df_summary, "variant_filter_summary.tsv")
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_SITE_FILTERS.rda")
    }
  }
)
