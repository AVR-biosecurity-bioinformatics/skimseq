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
      NULL
    )
    invisible(lapply(
      process_packages,
      library,
      character.only = TRUE,
      warn.conflicts = FALSE
    ))

    ### run code

    # Find files
    files <- list.files(pattern = "metrics.tsv.gz$", full.names = TRUE)

    # Define default binning rules for each metric. NA's will be estimated from files
    metric_specs <- tibble::tribble(
      ~RULE    , ~MIN , ~MAX , ~NBINS ,
      "QUAL"   ,    0 , NA   ,    100 ,
      "DP"     ,    0 , NA   ,    100 ,
      "ExcHet" ,    0 , 1    ,    100 ,
      "HWE"    ,    0 , 1    ,    100 ,
      "MAF"    ,    0 , 0.5  ,    100 ,
      "NS"     ,    0 , NA   ,    100 ,
      "CR"     ,    0 , 1    ,    100
    )

    # Function to subsample a few files and estimate max
    estimate_max_from_sample <- function(
      files,
      columns,
      sample_n = 5,
      seed = 1,
      buffer_frac = 0.05
    ) {
      set.seed(seed)
      sampled_files <- sample(files, min(sample_n, length(files)))

      max_list <- purrr::map(sampled_files, function(file) {
        df <- readr::read_tsv(
          file,
          show_col_types = FALSE,
          na = ".",
          col_select = dplyr::any_of(columns)
        )

        vapply(
          df,
          function(x) {
            x <- suppressWarnings(as.numeric(x))
            x <- x[is.finite(x)]
            if (length(x) == 0) NA_real_ else max(x)
          },
          numeric(1)
        )
      })

      max_df <- dplyr::bind_rows(lapply(max_list, as.list)) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::everything(),
            ~ if (all(is.na(.x))) {
              NA_real_
            } else {
              max(.x, na.rm = TRUE) * (1 + buffer_frac)
            }
          )
        )

      as.list(max_df)
    }

    cols_to_estimate <- metric_specs %>%
      dplyr::filter(is.na(MAX)) %>%
      dplyr::pull(RULE) %>%
      unique()

    sampled_max <- estimate_max_from_sample(
      files = files,
      columns = cols_to_estimate,
      sample_n = 5,
      buffer_frac = 0.10
    )

    metric_specs <- metric_specs %>%
      dplyr::mutate(
        MAX = ifelse(is.na(MAX), unlist(sampled_max[RULE]), MAX)
      )

    # Get just the per-pop metrics
    per_pop_metric_specs <- metric_specs %>%
      filter(RULE %in% c("NS", "MAF", "HWE", "ExcHet", "CR")) %>%
      dplyr::select(RULE, MIN, MAX, NBINS)

    # Function to read each chunked metrics file, and split into global and per_pop data frames
    split_metrics_one_chunk <- function(file) {
      df <- readr::read_tsv(file, show_col_types = FALSE, na = ".")

      # columns to extract from each summary file
      per_pop_prefixes <- c("NS", "MAF", "HWE", "ExcHet", "CR")
      per_pop_pattern <- paste0(
        "^(",
        paste(per_pop_prefixes, collapse = "|"),
        ")_"
      )

      per_pop_cols <- names(df)[stringr::str_detect(names(df), per_pop_pattern)]
      global_cols <- setdiff(
        names(df)[!names(df) %in% c("CHROM", "POS", "FILTER", "TYPE")],
        per_pop_cols
      )

      global_df <- df %>%
        select(CHROM, POS, FILTER, TYPE, all_of(global_cols)) %>%
        pivot_longer(
          cols = all_of(global_cols),
          names_to = "RULE",
          values_to = "VALUE"
        )

      per_pop_df <- df %>%
        select(CHROM, POS, FILTER, TYPE, all_of(per_pop_cols)) %>%
        pivot_longer(
          cols = all_of(per_pop_cols),
          names_to = "METRIC_POP",
          values_to = "VALUE"
        ) %>%
        tidyr::extract(
          METRIC_POP,
          into = c("RULE", "POP"),
          regex = "^(NS|MAF|HWE|ExcHet|CR)_(.+)$",
          remove = TRUE
        )

      list(global = global_df, per_pop = per_pop_df)
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

    # Function to summarise global and per_pop filter metrics for each chunk
    summarise_metric_chunk <- function(df, metric_specs, per_pop = FALSE) {
      fail_tag_fun <- if (per_pop) {
        function(type, rule) toupper(paste0("POP_", rule, "_FAIL"))
      } else {
        function(type, rule) toupper(paste0(rule, "_FAIL"))
      }

      df1 <- df %>%
        filter(!is.na(VALUE)) %>%
        left_join(metric_specs, by = "RULE") %>%
        mutate(
          ROW_ID = dplyr::row_number(),
          EXPECTED_FAIL = fail_tag_fun(TYPE, RULE)
        )

      # Expand FILTER tags once, then keep only exact matches to EXPECTED_FAIL
      fail_ids <- df1 %>%
        select(ROW_ID, FILTER, EXPECTED_FAIL) %>%
        tidyr::separate_longer_delim(FILTER, delim = ";") %>%
        mutate(
          FILTER = case_when(
            FILTER %in%
              c(
                "DP_LOWER_PERC_FAIL",
                "DP_MIN_FAIL",
                "DP_UPPER_PERC_FAIL"
              ) ~ "DP_FAIL",
            TRUE ~ FILTER
          )
        ) %>%
        filter(FILTER == EXPECTED_FAIL) %>%
        distinct(ROW_ID) %>%
        mutate(FILTER_CLASS = "FAIL")

      df2 <- df1 %>%
        select(-FILTER) %>%
        left_join(fail_ids, by = "ROW_ID") %>%
        mutate(
          FILTER_CLASS = dplyr::coalesce(FILTER_CLASS, "PASS")
        ) %>%
        group_by(RULE) %>%
        mutate(
          BIN = bin_values(VALUE, first(MIN), first(MAX), first(NBINS))
        ) %>%
        ungroup()

      if (per_pop) {
        df2 %>%
          count(RULE, POP, FILTER = FILTER_CLASS, TYPE, BIN, name = "COUNT")
      } else {
        df2 %>%
          count(RULE, FILTER = FILTER_CLASS, TYPE, BIN, name = "COUNT")
      }
    }

    # Process each chunk at a time
    global_list <- vector("list", length(files))
    per_pop_list <- vector("list", length(files))
    for (i in seq_along(files)) {
      #print(i)
      x <- split_metrics_one_chunk(files[[i]])
      if (nrow(x$global) > 0) {
        global_list[[i]] <- summarise_metric_chunk(
          x$global,
          metric_specs,
          per_pop = FALSE
        )
      }
      if (nrow(x$per_pop) > 0) {
        per_pop_list[[i]] <- summarise_metric_chunk(
          x$per_pop,
          per_pop_metric_specs,
          per_pop = TRUE
        )
      }
    }

    global_df <- bind_rows(global_list) %>%
      group_by(RULE, FILTER, TYPE, BIN) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop") %>%
      group_by(TYPE, RULE) %>%
      mutate(
        PROP = COUNT / sum(COUNT),
        n_grp = dplyr::n()
      ) %>%
      ungroup() %>%
      left_join(
        metric_specs %>% select(RULE, MIN, MAX),
        by = "RULE"
      )

    global_axis_df <- global_df %>%
      distinct(TYPE, RULE, MIN, MAX) %>%
      tidyr::pivot_longer(
        cols = c(MIN, MAX),
        names_to = "BOUND",
        values_to = "X"
      ) %>%
      mutate(PROP = 0, COUNT = 0)

    variant_types <- intersect(
      c("SNP", "INDEL", "REF", "ALL"),
      unique(global_df$TYPE)
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

      plot_df <- global_df %>%
        filter(TYPE == variant_type) %>%
        filter(is.finite(BIN), !is.na(BIN))

      axis_df <- global_axis_df %>%
        filter(TYPE == variant_type)

      if (nrow(plot_df) == 0) {
        next
      }

      global_qc_plots[[v]] <- ggplot(
        plot_df,
        aes(x = BIN, y = COUNT, fill = FILTER)
      ) +
        geom_blank(data = axis_df, aes(x = X, y = COUNT), inherit.aes = FALSE) +
        geom_col() +
        facet_wrap(~RULE, scales = "free") +
        scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
        scale_x_binned(n.breaks = 100, labels = thin_binned_labels(8)) +
        theme_classic() +
        labs(title = variant_type, x = NULL, y = "Number of sites") +
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
    per_pop_df <- bind_rows(per_pop_list) %>%
      group_by(RULE, POP, FILTER, TYPE, BIN) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop") %>%
      group_by(POP, TYPE, RULE) %>%
      mutate(PROP = COUNT / sum(COUNT)) %>%
      ungroup() %>%
      left_join(per_pop_metric_specs, by = "RULE")

    per_pop_axis_df <- per_pop_df %>%
      distinct(POP, RULE, TYPE, MIN, MAX) %>%
      tidyr::pivot_longer(
        cols = c(MIN, MAX),
        names_to = "BOUND",
        values_to = "X"
      ) %>%
      mutate(PROP = 0, COUNT = 0)

    #TODO: add back in  vlines for filter thresholds
    poprules <- expand_grid(unique(per_pop_df$RULE), unique(per_pop_df$TYPE))
    colnames(poprules) <- c("RULE", "TYPE")
    pop_qc_plots <- vector("list", length = nrow(poprules))

    for (p in 1:nrow(poprules)) {
      rule <- poprules$RULE[[p]]
      type <- poprules$TYPE[[p]]

      plot_df <- per_pop_df %>%
        filter(RULE == rule, TYPE == type) %>%
        filter(is.finite(BIN), !is.na(BIN))

      axis_df <- per_pop_axis_df %>%
        filter(RULE == rule)

      if (nrow(plot_df) == 0) {
        next
      }

      pop_qc_plots[[p]] <- ggplot(
        plot_df,
        aes(x = BIN, y = COUNT, fill = FILTER)
      ) +
        geom_blank(
          data = axis_df,
          aes(x = X, y = COUNT),
          inherit.aes = FALSE
        ) +
        geom_col() +
        facet_grid(POP ~ .) +
        scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
        scale_x_binned(
          n.breaks = 100,
          labels = thin_binned_labels(8)
        ) +
        theme_classic() +
        labs(
          title = paste(type, "per-population", rule),
          x = NULL,
          y = "Number of sites"
        ) +
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