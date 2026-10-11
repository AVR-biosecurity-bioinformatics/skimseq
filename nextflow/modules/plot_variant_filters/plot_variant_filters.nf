nextflow.enable.types = true

process PLOT_VARIANT_FILTERS {

    conda "${moduleDir}/environment.yml"

    input:
    histograms: Iterable<Path>

    stage:
    stageAs histograms, 'histograms/*'

    output:
    files("vcf_filters_*.pdf")

    script:
    '''
    #!/usr/bin/env bash
    set -euo pipefail

    unset R_LIBS R_LIBS_USER R_LIBS_SITE

    functions_r=$(command -v functions.R) || {
        echo "ERROR: functions.R not found on PATH" >&2
        exit 1
    }

    Rscript --vanilla - "$functions_r" <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)

    suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
    })

    source(args[[1L]])

    plot_counts <- read_variant_filter_histograms("histograms")

    # Population layout: two columns and up to three rows per page.
    populations_per_page <- 6L
    population_columns <- 2L

    # Preserve metric-specific PASS/FAIL classification.
    plot_counts[, FILTER := factor(
        FILTER,
        levels = c("PASS", "FAIL")
    )]

    setorder(plot_counts, RULE, POP, TYPE, BIN, FILTER)

    # Stack PASS and FAIL within each bin.
    plot_counts[
        ,
        `:=`(
            YMAX = cumsum(COUNT),
            YMIN = cumsum(COUNT) - COUNT
        ),
        by = .(RULE, POP, TYPE, BIN)
    ]

    rule_order <- c(
        "QUAL", "DP", "ExcHet", "HWE", "MAF", "NS", "CR"
    )

    plot_counts[, RULE := factor(
        RULE,
        levels = c(
            rule_order,
            sort(setdiff(
                as.character(unique(RULE)),
                rule_order
            ))
        )
    )]

    colours <- c(
        PASS = "#619CFF",
        FAIL = "#F8766D"
    )

    plot_theme <- theme_classic(base_size = 11) +
        theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(
                fill = "grey95",
                colour = "black"
            ),
            strip.text = element_text(size = 11),
            panel.spacing = grid::unit(1, "lines")
        )

    filter_plot <- function(
        x,
        variant_type,
        selected_rule,
        per_population = FALSE,
        page = NULL,
        pages = NULL,
        x_limits = NULL,
        y_limits = NULL
    ) {
        subtitle <- if (per_population) {
            paste0(
                "Population batch ", page, " of ", pages,
                "; shared axes across batches"
            )
        } else {
            "Stacked metric-specific PASS and FAIL counts"
        }

        p <- ggplot(
            x,
            aes(
                xmin = XMIN,
                xmax = XMAX,
                ymin = YMIN,
                ymax = YMAX,
                fill = FILTER
            )
        ) +
            geom_rect() +
            scale_fill_manual(
                values = colours,
                breaks = c("PASS", "FAIL"),
                drop = FALSE
            ) +
            scale_x_continuous(
                breaks = scales::breaks_pretty(n = 6),
                expand = expansion(mult = c(0, 0.02))
            ) +
            scale_y_continuous(
                breaks = scales::breaks_pretty(n = 4),
                labels = scales::label_number(
                    big.mark = ",",
                    scale_cut = scales::cut_short_scale()
                ),
                expand = expansion(mult = c(0, 0.05))
            ) +
            coord_cartesian(
                xlim = x_limits,
                ylim = y_limits
            ) +
            plot_theme +
            labs(
                title = paste(
                    variant_type,
                    if (per_population) {
                        "per-population filter diagnostics:"
                    } else {
                        "filter diagnostics:"
                    },
                    selected_rule
                ),
                subtitle = subtitle,
                caption = paste(
                    "PASS/FAIL refers to this metric only,",
                    "not the overall site FILTER status."
                ),
                x = selected_rule,
                y = "Number of sites"
            )

        if (per_population) {
            p + facet_wrap(
                ~POP,
                ncol = population_columns,
                scales = "fixed"
            )
        } else {
            p
        }
    }

    ordered_types <- function(x) {
        preferred <- c("SNP", "INDEL", "REF", "ALL")
        observed <- unique(as.character(x))

        c(
            intersect(preferred, observed),
            sort(setdiff(observed, preferred))
        )
    }

    empty_page <- function(message) {
        plot.new()
        text(0.5, 0.5, message)
    }

    variant_types <- ordered_types(plot_counts$TYPE)

    for (variant_type in variant_types) {
        file_type <- tolower(variant_type)

        # Global: separate file per variant type,
        # with one statistic per page.
        pdf(
            paste0("vcf_filters_", file_type, "_global.pdf"),
            width = 11,
            height = 8
        )

        present_rules <- plot_counts[
            POP == "." & TYPE == variant_type,
            unique(as.character(RULE))
        ]

        rules <- intersect(
            levels(plot_counts$RULE),
            present_rules
        )

        if (length(rules)) {
            for (selected_rule in rules) {
                x <- plot_counts[
                    POP == "." &
                    TYPE == variant_type &
                    as.character(RULE) == selected_rule
                ]

                print(filter_plot(
                    x,
                    variant_type,
                    selected_rule
                ))

                rm(x)
                invisible(gc())
            }
        } else {
            empty_page(paste(
                "No global histogram data for",
                variant_type
            ))
        }

        dev.off()

        # Population: separate file per variant type.
        # Each statistic spans as many six-population pages as needed.
        pdf(
            paste0("vcf_filters_", file_type, "_perpop.pdf"),
            width = 11,
            height = 8
        )

        present_rules <- plot_counts[
            POP != "." & TYPE == variant_type,
            unique(as.character(RULE))
        ]

        rules <- intersect(
            levels(plot_counts$RULE),
            present_rules
        )

        if (length(rules)) {
            for (selected_rule in rules) {
                x <- plot_counts[
                    POP != "." &
                    TYPE == variant_type &
                    as.character(RULE) == selected_rule
                ]

                populations <- sort(unique(as.character(x$POP)))

                batches <- split(
                    populations,
                    ceiling(
                        seq_along(populations) / populations_per_page
                    )
                )

                # Keep scales identical across all batches of this metric.
                x_limits <- range(c(x$XMIN, x$XMAX))
                y_limits <- c(0, max(x$YMAX))

                if (y_limits[[2L]] == 0) {
                    y_limits[[2L]] <- 1
                }

                for (page in seq_along(batches)) {
                    selected_populations <- batches[[page]]

                    page_df <- copy(
                        x[POP %in% selected_populations]
                    )

                    page_df[, POP := factor(
                        POP,
                        levels = selected_populations
                    )]

                    print(filter_plot(
                        page_df,
                        variant_type,
                        selected_rule,
                        per_population = TRUE,
                        page = page,
                        pages = length(batches),
                        x_limits = x_limits,
                        y_limits = y_limits
                    ))

                    rm(page_df)
                    invisible(gc())
                }

                rm(x)
                invisible(gc())
            }
        } else {
            empty_page(paste(
                "No per-population histogram data for",
                variant_type
            ))
        }

        dev.off()
    }
    RSCRIPT
    '''
}