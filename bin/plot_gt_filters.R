#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]

    #snp_filtering parameters
    params.gt_qual <- args[3]
    params.gt_dp_min <- args[4]
    params.gt_dp_max <- args[5]

    # TESTING
    #params.gt_qual <- 20
    #params.gt_dp_min <- 1
    #params.gt_dp_max <- 100

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

    # Create filtering parameter table
    # fmt: skip
    parameter_table <- tibble::tribble(
      ~filter, ~type, ~lower, ~upper,
      "GQ", "SNP", params.gt_qual, NA,
      "DP", "SNP", params.gt_dp_min, params.gt_dp_max,
      "GQ", "INDEL", params.gt_qual, NA,
      "DP", "INDEL", params.gt_dp_min, params.gt_dp_max,
      "GQ", "NO_VARIATION", params.gt_qual, NA,
      "DP", "NO_VARIATION", params.gt_dp_min, params.gt_dp_max
    ) %>%
      dplyr::mutate(
        upper = as.numeric(upper),
        lower = as.numeric(lower)
      )

    # List table files
    table_files <- list.files(pattern = "\\.table.gz$")

    # loop through parameter table
    plot_list <- vector("list", length = nrow(parameter_table))
    for (i in 1:nrow(parameter_table)) {
      filter_to_select <- parameter_table$filter[i]
      type_to_select <- parameter_table$type[i]

      df <- read_tsv(
        table_files,
        col_select = c("TYPE", "FILTER", "COUNTS", filter_to_select)
      ) %>%
        dplyr::filter(TYPE == type_to_select) %>%
        dplyr::filter(!FILTER == "MISSING") %>%
        dplyr::mutate(
          label = ifelse(
            stringr::str_detect(FILTER, filter_to_select),
            "FAIL",
            "PASS"
          )
        )

      plot_list[[i]] <- df %>%
        ggplot(aes(x = !!sym(filter_to_select), y = COUNTS, fill = label)) +
        geom_col() +
        #geom_density()+
        geom_vline(
          data = parameter_table %>% dplyr::slice(i),
          aes(xintercept = lower),
          lty = "dashed"
        ) +
        geom_vline(
          data = parameter_table %>% dplyr::slice(i),
          aes(xintercept = upper),
          lty = "dashed"
        ) +
        scale_fill_manual(values = c("PASS" = "#619CFF", "FAIL" = "#F8766D")) +
        labs(
          x = paste0(type_to_select, " ", filter_to_select),
          y = "Count"
        ) +
        theme_classic() +
        theme(legend.position = "none")
    }

    gg.genotype_qc_dist <- wrap_plots(plot_list)

    # Write out plots
    pdf("genotype_filtering_qc.pdf", width = 11, height = 8)
    plot(gg.genotype_qc_dist)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_VARIANT_QC.rda")
    }
  }
)
