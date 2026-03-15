#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]
    outname <- args[3]

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

    # columns to extract from each summary file
    global_cols <- c(
      "CHROM",
      "POS",
      "FILTER",
      "QUAL",
      "DP",
      "CR",
      "NS",
      "MAF",
      "HWE",
      "TYPE",
      "ExcHet",
      "DIST_INDEL"
    )
    per_pop_prefixes <- c("NS", "MAF", "HWE", "ExcHet")
    per_pop_pattern <- paste0(
      "^(",
      paste(per_pop_prefixes, collapse = "|"),
      ")_"
    )

    # Function to read each chunked metrics file, and split into global and per_pop data frames
    split_metrics_one_chunk <- function(file) {
      df <- readr::read_tsv(file, show_col_types = FALSE, na = ".") %>%
        mutate(
          VARIANT_TYPE = dplyr::case_when(
            TYPE == "SNP" ~ "snp",
            TYPE == "INDEL" ~ "indel",
            TYPE == "REF" ~ "invariant",
            TRUE ~ tolower(TYPE)
          )
        )

      per_pop_cols <- names(df)[stringr::str_detect(names(df), per_pop_pattern)]
      global_cols <- setdiff(names(df), per_pop_cols)

      global_df <- df %>%
        select(all_of(global_cols))

      per_pop_df <- df %>%
        select(CHROM, POS, FILTER, TYPE, VARIANT_TYPE, all_of(per_pop_cols)) %>%
        pivot_longer(
          cols = all_of(per_pop_cols),
          names_to = "METRIC_POP",
          values_to = "VALUE"
        ) %>%
        tidyr::extract(
          METRIC_POP,
          into = c("RULE", "POP"),
          regex = "^(NS|MAF|HWE|ExcHet)_(.+)$",
          remove = TRUE
        )

      list(global = global_df, per_pop = per_pop_df)
    }

    # Define binning rules for each metric
    # TODO  - estimate these histogram bin sizes from one of the chunk files
    metric_specs <- tibble::tribble(
      ~RULE        , ~COLUMN      , ~FAIL_PATTERN                                 , ~MIN , ~MAX    , ~NBINS ,
      "QUAL"       , "QUAL"       , "SNP_QUAL_FAIL|INDEL_QUAL_FAIL|INV_QUAL_FAIL" ,    0 ,  5000   ,    100 ,
      "DP"         , "DP"         , "SNP_DP_FAIL|INDEL_DP_FAIL|INV_DP_FAIL"       ,    0 , 50000   ,    100 ,
      "ExcHet"     , "ExcHet"     , "SNP_EH_FAIL|INDEL_EH_FAIL"                   ,    0 ,     1   ,    100 ,
      "HWE"        , "HWE"        , "SNP_HWE_FAIL|INDEL_HWE_FAIL"                 ,    0 ,     1   ,    100 ,
      "MAF"        , "MAF"        , "SNP_MAF_FAIL|INDEL_MAF_FAIL"                 ,    0 ,     0.5 ,    100 ,
      "NS"         , "NS"         , "SNP_NS_FAIL|INDEL_NS_FAIL|INV_NS_FAIL"       ,    0 ,  2000   ,    100 ,
      "CR"         , "CR"         , "SNP_CR_FAIL|INDEL_CR_FAIL|INV_CR_FAIL"       ,    0 ,     1   ,    100 ,
      "DIST_INDEL" , "DIST_INDEL" , "SNP_DIST_INDEL_FAIL"                         ,    0 ,   100   ,    100
    )
    per_pop_metric_specs <- tibble::tribble(
      ~RULE    , ~MIN , ~MAX   , ~NBINS ,
      "NS"     ,    0 , 2000   ,    100 ,
      "MAF"    ,    0 ,    0.5 ,    100 ,
      "HWE"    ,    0 ,    1   ,    100 ,
      "ExcHet" ,    0 ,    1   ,    100
    )

    # Function to bin statistic values to reduce file size
    bin_values <- function(x, min_val, max_val, nbins) {
      x <- as.numeric(x)
      x[x < min_val] <- min_val
      x[x > max_val] <- max_val
      breaks <- seq(min_val, max_val, length.out = nbins + 1)
      mids <- head(breaks, -1) + diff(breaks) / 2
      idx <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      mids[idx]
    }

    # Function to summarise metrics for each chunk
    summarise_global_chunk <- function(global_df, metric_specs) {
      purrr::map_dfr(seq_len(nrow(metric_specs)), function(i) {
        spec <- metric_specs[i, ]

        vals <- global_df[[spec$COLUMN]]
        keep <- !is.na(vals)

        if (!any(keep)) {
          return(NULL)
        }

        tibble(
          RULE = spec$RULE,
          VARIANT_TYPE = global_df$VARIANT_TYPE[keep],
          FILTER = if_else(
            stringr::str_detect(global_df$FILTER[keep], spec$FAIL_PATTERN),
            "FAIL",
            "PASS"
          ),
          BIN = bin_values(vals[keep], spec$MIN, spec$MAX, spec$NBINS)
        ) %>%
          count(RULE, FILTER, VARIANT_TYPE, BIN, name = "COUNT")
      })
    }

    summarise_per_pop_chunk <- function(per_pop_df, per_pop_metric_specs) {
      per_pop_df %>%
        left_join(per_pop_metric_specs, by = "RULE") %>%
        filter(!is.na(VALUE)) %>%
        mutate(
          FILTER = if_else(stringr::str_detect(FILTER, "POP"), "FAIL", "PASS"),
          BIN = purrr::pmap_dbl(
            list(VALUE, MIN, MAX, NBINS),
            ~ bin_values(..1, ..2, ..3, ..4)
          )
        ) %>%
        count(RULE, POP, FILTER, VARIANT_TYPE, BIN, name = "COUNT")
    }

    files <- list.files(pattern = "metrics.tsv.gz$", full.names = TRUE)

    # Process each chunk at a time
    global_list <- vector("list", length(files))
    per_pop_list <- vector("list", length(files))

    for (i in seq_along(files)) {
      x <- split_metrics_one_chunk(files[[i]])

      global_list[[i]] <- summarise_global_chunk(x$global, metric_specs)
      per_pop_list[[i]] <- summarise_per_pop_chunk(
        x$per_pop,
        per_pop_metric_specs
      )

      rm(x)
      gc(FALSE)
    }

    global_df <- bind_rows(global_list) %>%
      group_by(RULE, FILTER, VARIANT_TYPE, BIN) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop") %>%
      group_by(VARIANT_TYPE, RULE) %>%
      mutate(
        PROP = COUNT / sum(COUNT),
        n_grp = dplyr::n()
      ) %>%
      ungroup()

    per_pop_df <- bind_rows(per_pop_list) %>%
      group_by(RULE, POP, FILTER, VARIANT_TYPE, BIN) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop") %>%
      group_by(POP, VARIANT_TYPE, RULE) %>%
      mutate(PROP = COUNT / sum(COUNT)) %>%
      ungroup()

    variant_types <- factor(
      unique(global_df$VARIANT_TYPE),
      levels = c("snp", "indel", "invariant", "all")
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

    # Create global qc plots
    #TODO: add back in  vlines for filter thresholds
    global_qc_plots <- vector("list", length = length(variant_types))
    for (v in 1:length(variant_types)) {
      variant_type <- variant_types[v]
      global_qc_plots[[v]] <- global_df_aug %>%
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
    pdf(paste0(outname, "_global.pdf"), width = 11, height = 8)
    purrr::walk(global_qc_plots, plot)
    try(dev.off(), silent = TRUE)

    # Create per-pop qc plots
    #TODO: add back in  vlines for filter thresholds
    poprules <- unique(per_pop_df$RULE)
    pop_qc_plots <- vector("list", length = length(poprules))

    for (p in 1:length(poprules)) {
      rule <- poprules[[p]]
      pop_qc_plots[[p]] <- per_pop_df %>%
        filter(RULE == rule) %>%
        ggplot(aes(x = BIN, y = PROP, fill = FILTER)) +
        geom_col() +
        facet_grid(POP ~ paste0(VARIANT_TYPE, ":", RULE), scales = "free") +
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
    pdf(paste0(outname, "_perpop.pdf"), width = 11, height = 8)
    purrr::walk(pop_qc_plots, plot)
    try(dev.off(), silent = TRUE)

    # Create a joint table of the summary files
    # List filtering summary files

    # TODO: this needs to be summarised by variant type
    #summary_files <- list.files(pattern = "filter_summary.tsv$")
    #summary_files <- summary_files[file.size(summary_files) > 0]
    #df_summary <- readr::read_tsv(
    #  summary_files,
    #  col_names = c("FILTER", "COUNT"),
    #  col_types = c("cn")
    #) %>%
    #  group_by(FILTER) %>%
    #  summarise(COUNT = sum(COUNT))

    # Write out summary file
    #write_tsv(df_summary, paste0(outname, ".tsv"))
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_VCF_FILTERS.rda")
    }
  }
)
