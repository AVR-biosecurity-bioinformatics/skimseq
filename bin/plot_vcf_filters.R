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

    # Function to Normalise variant names
    normalise_type <- function(x) {
      dplyr::case_when(
        x == "SNP" ~ "snp",
        x == "INDEL" ~ "indel",
        x == "REF" ~ "invariant",
        TRUE ~ tolower(x)
      )
    }

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
    summarise_one_chunk <- function(file, metric_specs) {
      df <- readr::read_tsv(file, show_col_types = FALSE) %>%
        mutate(VARIANT_TYPE = normalise_type(TYPE))

      purrr::map_dfr(seq_len(nrow(metric_specs)), function(i) {
        spec <- metric_specs[i, ]

        vals <- df[[spec$COLUMN]]
        keep <- !is.na(vals)

        if (!any(keep)) {
          return(NULL)
        }

        tibble(
          RULE = spec$RULE,
          VARIANT_TYPE = df$VARIANT_TYPE[keep],
          FILTER = if_else(
            stringr::str_detect(df$FILTER[keep], spec$FAIL_PATTERN),
            "FAIL",
            "PASS"
          ),
          BIN = bin_values(vals[keep], spec$MIN, spec$MAX, spec$NBINS)
        ) %>%
          count(RULE, FILTER, VARIANT_TYPE, BIN, name = "COUNT")
      })
    }

    files <- list.files(pattern = "tagged_metrics.tsv.gz$", full.names = TRUE)

    chunk_summaries <- lapply(
      files,
      summarise_one_chunk,
      metric_specs = metric_specs
    )

    # Read in all tables and combine
    df <- bind_rows(chunk_summaries) %>%
      group_by(RULE, FILTER, VARIANT_TYPE, BIN) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop") %>%
      group_by(VARIANT_TYPE, RULE) %>%
      mutate(PROP = COUNT / sum(COUNT)) %>%
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
    pdf(paste0(outname, ".pdf"), width = 11, height = 8)
    purrr::walk(variant_qc_plots, plot)
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
