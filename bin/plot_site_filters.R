#!/usr/bin/env Rscript

tryCatch(
  {
    args <- commandArgs(trailingOnly = TRUE)

    projectDir <- args[1]
    params.rdata <- args[2]

    #snp_filtering parameters
    params.snp_qd <- args[3]
    params.snp_qual <- args[4]
    params.snp_sor <- args[5]
    params.snp_fs <- args[6]
    params.snp_mq <- args[7]
    params.snp_mqrs <- args[8]
    params.snp_rprs <- args[9]
    params.snp_maf <- args[10]
    params.snp_mac <- args[11]
    params.snp_eh <- args[12]
    params.snp_dp_min <- args[13]
    params.snp_dp_max <- args[14]

    #indel_filtering parameters
    params.indel_qd <- args[15]
    params.indel_qual <- args[16]
    params.indel_fs <- args[17]
    params.indel_rprs <- args[18]
    params.indel_maf <- args[19]
    params.indel_mac <- args[20]
    params.indel_eh <- args[21]
    params.indel_dp_min <- args[22]
    params.indel_dp_max <- args[23]

    #invariant_filtering parameters
    params.inv_dp_min <- args[24]
    params.inv_dp_max <- args[25]

    # Max missing
    params.snp_max_missing <- args[26]
    params.indel_max_missing <- args[27]
    params.inv_max_missing <- args[28]

    # TESTING
    #params.snp_qd <- 2.0
    #params.snp_qual <- 30.0
    #params.snp_sor <- 3.0
    #params.snp_fs <- 60.0
    #params.snp_mq <- 40.0
    #params.snp_mqrs <- -12.5
    #params.snp_rprs <- -8.0
    #params.snp_maf <- 0.05
    #params.snp_eh <- 54.69
    #params.snp_dp_min <- 6
    #params.snp_dp_max <- 1500
    #params.indel_qd <- 2.0
    #params.indel_qual <- 30.0
    #params.indel_fs <- 200.0
    #params.indel_rprs <- -20.0
    #params.indel_maf <- 0.05
    #params.indel_eh <- 54.69
    #params.indel_dp_min <- 6
    #params.indel_dp_max <- 1500
    #params.inv_dp_min <- 6
    #params.inv_dp_max <- 1500
    #params.max_missing <- 0.05

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
      "QD", "snp", params.snp_qd, NA,
      "QUAL", "snp", params.snp_qual, NA,
      "SOR", "snp", params.snp_sor, NA,
      "FS", "snp", params.snp_fs, NA,
      "MQ", "snp", params.snp_mq, NA,
      "MQRankSum", "snp", params.snp_mqrs, NA,
      "ReadPosRankSum", "snp", params.snp_rprs, NA,
      "MAF", "snp", params.snp_maf, NA,
      "MAC", "snp", params.snp_mac, NA,
      "ExcessHet", "snp", params.snp_eh, NA,
      "DP", "snp", params.snp_dp_min, params.snp_dp_max,
      "QD", "indel", params.indel_qd, NA,
      "QUAL", "indel", params.indel_qual, NA,
      "FS", "indel", params.indel_fs, NA,
      "ReadPosRankSum", "indel", params.indel_rprs, NA,
      "MAF", "indel", params.indel_maf, NA,
      "MAC", "indel", params.indel_mac, NA,
      "ExcessHet", "indel", params.indel_eh, NA,
      "DP", "indel", params.indel_dp_min, params.indel_dp_max,
      "DP", "invariant", params.inv_dp_min, params.inv_dp_max,
      "F_MISSING", "snp", params.snp_max_missing, NA,
      "F_MISSING", "indel", params.indel_max_missing, NA,
      "F_MISSING", "invariant", params.inv_max_missing, NA,
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

      table_to_read <- table_files[stringr::str_detect(
        table_files,
        type_to_select
      )]

      df <- read_tsv(
        table_to_read,
        col_select = c("FILTER", filter_to_select),
        col_types = c("cn")
      ) %>%
        dplyr::mutate(
          label = ifelse(
            stringr::str_detect(FILTER, filter_to_select),
            "FAIL",
            "PASS"
          )
        )

      plot_list[[i]] <- df %>%
        ggplot(aes(x = !!sym(filter_to_select), fill = label)) +
        geom_histogram() +
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

    gg.snp_qc_dist <- wrap_plots(plot_list)

    # Write out plots
    pdf("variant_filtering_qc.pdf", width = 11, height = 8)
    plot(gg.snp_qc_dist)
    try(dev.off(), silent = TRUE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "PLOT_VARIANT_QC.rda")
    }
  }
)
