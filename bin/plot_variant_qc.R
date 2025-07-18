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
    params.snp_eh <- args[11]
    params.snp_dp_min <- args[12]
    params.snp_dp_max <- args[13]

    #indel_filtering parameters
    params.indel_qd <- args[14]
    params.indel_qual <- args[15]
    params.indel_fs <- args[16]
    params.indel_rprs <- args[17]
    params.indel_maf <- args[18]
    params.indel_eh <- args[19]
    params.indel_dp_min <- args[20]
    params.indel_dp_max <- args[21]

    #invariant_filtering parameters
    params.inv_dp_min <- args[22]
    params.inv_dp_max <- args[23]

    # Max missing
    params.max_missing <- args[24]

    # TESTING
    #snp_filtering parameters
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

    #indel_filtering parameters
    #params.indel_qd <- 2.0
    #params.indel_qual <- 30.0
    #params.indel_fs <- 200.0
    #params.indel_rprs <- -20.0
    #params.indel_maf <- 0.05
    #params.indel_eh <- 54.69
    #params.indel_dp_min <- 6
    #params.indel_dp_max <- 1500

    #invariant_filtering parameters
    #params.inv_dp_min <- 6
    #params.inv_dp_max <- 1500

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
    # TODO - where doe the LowQual filter come from???

    # List table files
    table_files <- list.files(pattern = "\\.table.gz$")

    # read in all table files and turn into long format

    # TODO: Do this in bash
    snp_qc <- readr::read_tsv(table_files) %>%
      dplyr::select(-NS) %>% # NS is the inverse of F_Missing
      tidyr::pivot_longer(-c("CHROM", "POS", "TYPE", "FILTER")) %>%
      tidyr::separate_longer_delim(cols = FILTER, delim = ",")

    # Update names of filters to match names of parameters
    trans_table <- snp_qc %>%
      dplyr::select(FILTER) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        FILTER_NEW = case_when(
          FILTER %in% c("DPmin", "DPmax") ~ "DP",
          TRUE ~ FILTER
        )
      )

    # Create filtering parameter table
    # fmt: skip
    parameter_table <- tibble::tribble(
      ~name, ~TYPE, ~lower, ~upper,
      "QD", "SNP", params.snp_qd, NA,
      "QUAL", "SNP", params.snp_qual, NA,
      "SOR", "SNP", params.snp_sor, NA,
      "FS", "SNP", params.snp_fs, NA,
      "MQ", "SNP", params.snp_mq, NA,
      "MQRankSum", "SNP", params.snp_mqrs, NA,
      "ReadPosRankSum", "SNP", params.snp_rprs, NA,
      "AF", "SNP", params.snp_maf, NA,
      "ExcessHet", "SNP", params.snp_eh, NA,
      "DP", "SNP", params.snp_dp_min, params.snp_dp_max,
      "QD", "INDEL", params.indel_qd, NA,
      "QUAL", "INDEL", params.indel_qual, NA,
      "FS", "INDEL", params.indel_fs, NA,
      "ReadPosRankSum", "INDEL", params.indel_rprs, NA,
      "AF", "INDEL", params.indel_maf, NA,
      "ExcessHet", "INDEL", params.indel_eh, NA,
      "DP", "INDEL", params.indel_dp_min, params.indel_dp_max,
      "DP", "NO_VARIATION", params.inv_dp_min, params.inv_dp_max,
      "F_MISSING", "SNP", params.max_missing, NA,
      "F_MISSING", "INDEL", params.max_missing, NA,
      "F_MISSING", "NO_VARIATION", params.max_missing, NA,
    ) %>%
      dplyr::mutate(
        upper = as.numeric(upper),
        lower = as.numeric(lower)
      )

    # Do some reformatting to make sure SNPS arent counted twice
    snp_qc_dist <- snp_qc %>%
      dplyr::left_join(trans_table) %>%
      dplyr::mutate(
        pass = case_when(
          FILTER_NEW == name & !FILTER_NEW == "PASS" ~ 0,
          TRUE ~ 1
        )
      ) %>%
      dplyr::group_by(CHROM, POS, TYPE) %>%
      dplyr::slice_min(pass) %>%
      dplyr::mutate(
        label = case_when(
          pass == 1 ~ "PASS",
          pass == 0 ~ "FAIL"
        )
      )

    gg.snp_qc_dist <- snp_qc_dist %>%
      dplyr::mutate(
        TYPE = factor(
          TYPE,
          levels = c(
            "SNP",
            "INDEL",
            "NO_VARIATION"
          )
        )
      ) %>%
      ggplot(aes(x = value, fill = label)) +
      geom_histogram() +
      geom_vline(
        data = parameter_table,
        aes(xintercept = lower),
        lty = "dashed"
      ) +
      geom_vline(
        data = parameter_table,
        aes(xintercept = upper),
        lty = "dashed"
      ) +
      facet_grid(TYPE ~ name, scales = "free")

    # Write out plots
    pdf("variant_filtering_qc.pdf", width = 11, height = 8)
    plot(gg.snp_qc_dist)
    try(dev.off(), silent = TRUE)

    # TODO: create upset plot

    # TODO: Write out summary files
    readr::read_tsv(table_files) %>%
      dplyr::arrange(CHROM, POS) %>%
      dplyr::group_by(TYPE, FILTER) %>%
      dplyr::summarise(n = n_distinct(CHROM, POS)) %>%
      readr::write_tsv("variant_filtering_summary.tsv", col_names = FALSE)
  },
  finally = {
    ### save R environment if script throws error code
    if (params.rdata == "true") {
      save.image(file = "CREATE_INTERVALS.rda")
    }
  }
)
