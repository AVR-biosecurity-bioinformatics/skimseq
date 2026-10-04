process PLOT_VARIANT_FILTERS {

    conda "${moduleDir}/environment.yml"

    input:
    path 'histograms/*'

    output:
    path('vcf_filters_*.pdf'), emit: plots

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
    source(args[[1L]])

    suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
    })

    plot_counts <- read_variant_filter_histograms("histograms")

    # Stack PASS and FAIL counts within each bin.
    plot_counts[
        ,
        FILTER := factor(FILTER, levels = c("PASS", "FAIL"))
    ]

    setorder(plot_counts, RULE, POP, TYPE, BIN, FILTER)

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

    plot_counts[
        ,
        RULE := factor(
            RULE,
            levels = c(
                rule_order,
                setdiff(as.character(unique(RULE)), rule_order)
            )
        )
    ]

    colours <- c(
        PASS = "#619CFF",
        FAIL = "#F8766D"
    )

    plot_theme <- theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    # Global: one page per variant type, faceted by metric.
    global_df <- plot_counts[POP == "."]

    pdf("vcf_filters_global.pdf", width = 11, height = 8)

    if (nrow(global_df)) {
        preferred_types <- c("SNP", "INDEL", "REF", "ALL")

        types <- c(
            intersect(preferred_types, unique(global_df$TYPE)),
            setdiff(unique(global_df$TYPE), preferred_types)
        )

        for (variant_type in types) {
            x <- global_df[TYPE == variant_type]

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
                facet_wrap(~RULE, scales = "free") +
                scale_fill_manual(values = colours) +
                scale_x_continuous(
                    breaks = scales::breaks_pretty(n = 8)
                ) +
                plot_theme +
                labs(
                    title = variant_type,
                    x = NULL,
                    y = "Number of sites"
                )

            print(p)
        }
    } else {
        plot.new()
        text(0.5, 0.5, "No global histogram data")
    }

    dev.off()

    # Per-population: one page per metric and variant type.
    per_pop_df <- plot_counts[POP != "."]

    pdf("vcf_filters_perpop.pdf", width = 11, height = 8)

    if (nrow(per_pop_df)) {
        poprules <- unique(per_pop_df[, .(RULE, TYPE)])

        for (i in seq_len(nrow(poprules))) {
            selected_rule <- as.character(poprules$RULE[[i]])
            variant_type <- poprules$TYPE[[i]]

            x <- per_pop_df[
                as.character(RULE) == selected_rule &
                TYPE == variant_type
            ]

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
                facet_grid(POP ~ .) +
                scale_fill_manual(values = colours) +
                scale_x_continuous(
                    breaks = scales::breaks_pretty(n = 8)
                ) +
                plot_theme +
                labs(
                    title = paste(
                        variant_type,
                        "per-population",
                        selected_rule
                    ),
                    x = NULL,
                    y = "Number of sites"
                )

            print(p)
        }
    } else {
        plot.new()
        text(0.5, 0.5, "No per-population histogram data")
    }

    dev.off()
    RSCRIPT
    '''
}