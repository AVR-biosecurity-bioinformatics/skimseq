process PLOT_SAMPLE_FILTERS {
    tag "${missing_summary}"
    conda "${moduleDir}/environment.yml"

    input:
    path(missing_summary)

    output: 
    path("*.pdf"),               emit: plots
    path("sample_missing.tsv"),  emit: sample_missing_tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    Rscript --vanilla - <<'RSCRIPT'
        suppressPackageStartupMessages({
            library(ggplot2)
            library(readr)
        })

        threshold <- as.numeric(
            "${params.vcf_sample_max_missing}"
        )

        if (
            length(threshold) != 1L ||
            !is.finite(threshold) ||
            threshold < 0 ||
            threshold > 1
        ) {
            stop(
                "vcf_sample_max_missing must be a single value between 0 and 1; got: ",
                "${params.vcf_sample_max_missing}"
            )
        }

        sample_metrics <- readr::read_tsv(
            "${missing_summary}",
            show_col_types = FALSE,
            progress = FALSE
        )

        required_columns <- c(
            "MISSING_FRACTION"
        )

        missing_columns <- setdiff(
            required_columns,
            names(sample_metrics)
        )

        if (length(missing_columns) > 0L) {
            stop(
                "Missing required column(s) in ${missing_summary}: ",
                paste(missing_columns, collapse = ", ")
            )
        }

        sample_metrics\$MISSING_FRACTION <- suppressWarnings(
            as.numeric(sample_metrics\$MISSING_FRACTION)
        )

        invalid_missingness <- (
            !is.finite(sample_metrics\$MISSING_FRACTION) |
            sample_metrics\$MISSING_FRACTION < 0 |
            sample_metrics\$MISSING_FRACTION > 1
        )

        if (any(invalid_missingness)) {
            warning(
                "Ignoring ",
                sum(invalid_missingness),
                " sample(s) with invalid missing-data fractions"
            )
        }

        sample_metrics\$FILTER <- ifelse(
            invalid_missingness,
            NA_character_,
            ifelse(
                sample_metrics\$MISSING_FRACTION > threshold,
                "FAIL",
                "PASS"
            )
        )

        sample_metrics\$FILTER <- factor(
            sample_metrics\$FILTER,
            levels = c("PASS", "FAIL")
        )

        sample_plot <- ggplot(sample_metrics[!invalid_missingness, , drop = FALSE ],
            aes(x = MISSING_FRACTION, fill = FILTER )) +
            geom_histogram(
                bins = 30,
                colour = "white",
                linewidth = 0.2
            ) +
            geom_vline(
                xintercept = threshold,
                linetype = "dashed"
            ) +
            scale_fill_manual(
                values = c(
                    PASS = "#619CFF",
                    FAIL = "#F8766D"
                ),
                drop = FALSE
            ) +
            labs(
                x = "Missing-data fraction",
                y = "Number of samples"
            ) +
            theme_classic() +
            theme(
                legend.position = "none"
            )

        readr::write_tsv(
            sample_metrics,
            "sample_missing.tsv"
        )

        ggsave(
            filename = "sample_filtering_qc.pdf",
            plot = sample_plot,
            device = "pdf",
            width = 11,
            height = 8,
            units = "in"
        )
    RSCRIPT
    """
}