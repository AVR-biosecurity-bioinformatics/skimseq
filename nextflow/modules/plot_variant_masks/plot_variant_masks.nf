nextflow.enable.types = true

process PLOT_VARIANT_MASKS {

    conda "${moduleDir}/environment.yml"

    input:
    histograms: Iterable<Path>

    stage:
    stageAs histograms, 'histograms/*'

    output:
    files("vcf_masks_*.pdf")

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

    # Reading, validation and aggregation are handled by functions.R.
    mask_data <- read_variant_mask_histograms("histograms")

    plot_steps <- mask_data$steps
    record_totals <- mask_data$records

    rm(mask_data)
    invisible(gc())

    # Keep the union for count plots, but not distribution panels.
    plot_steps <- plot_steps[MASK != "AnyMask"]

    rule_order <- c(
        "QUAL", "DP", "ExcHet", "HWE", "MAF", "NS", "CR"
    )

    plot_steps[, RULE := factor(
        RULE,
        levels = c(
            rule_order,
            sort(setdiff(unique(RULE), rule_order))
        )
    )]

    present_masks <- unique(as.character(plot_steps$MASK))

    mask_order <- c(
        intersect("Unmasked", present_masks),
        sort(setdiff(present_masks, "Unmasked"))
    )

    # Readable labels for hash-free FILTER IDs.
    mask_labels <- setNames(
        sub("_FAIL$", "", sub("^MASK_", "", mask_order)),
        mask_order
    )

    colours <- setNames(
        grDevices::hcl.colors(
            max(1L, length(mask_order)),
            palette = "Dark 3"
        )[seq_along(mask_order)],
        mask_order
    )

    if ("Unmasked" %in% names(colours)) {
        colours["Unmasked"] <- "#619CFF"
    }

    plot_steps[, MASK := factor(MASK, levels = mask_order)]
    setorder(plot_steps, RULE, POP, TYPE, MASK, POSITION)

    plot_theme <- theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(
                fill = "grey95",
                colour = "black"
            )
        )

        distribution_plot <- function(
                x,
                variant_type,
                selected_rule,
                per_population = FALSE
            ) {
                # Work on a copy to avoid modifying the caller's data.table.
                x <- copy(x)

                page_masks <- mask_order[
                    mask_order %in% as.character(unique(x$MASK))
                ]

                x[, MASK := factor(
                    as.character(MASK),
                    levels = page_masks
                )]

                if (per_population) {
                    x[, POP := factor(
                        as.character(POP),
                        levels = sort(unique(as.character(POP)))
                    )]
                }

                # Complete bins using only positions belonging to the current rule.
                # Using plot_steps$POSITION here would mix the scales of unrelated
                # metrics such as DP, QUAL, HWE and MAF.
                rule_positions <- sort(unique(x$POSITION))

                x <- x[
                    CJ(
                        MASK = unique(x$MASK),
                        POP = unique(x$POP),
                        POSITION = rule_positions,
                        unique = TRUE
                    ),
                    on = .(MASK, POP, POSITION)
                ]

                x[is.na(PERCENT), PERCENT := 0]
                setorder(x, MASK, POP, POSITION)

                duplicate_bins <- x[
                    ,
                    .N,
                    by = .(MASK, POP, POSITION)
                ][N > 1L]

                if (nrow(duplicate_bins)) {
                    stop(
                        "Repeated POSITION values within MASK/POP groups. ",
                        "Check whether functions.R returns histogram bins or ",
                        "expanded step coordinates before using geom_line()."
                    )
                }

                p <- ggplot(
                    x,
                    aes(
                        x = POSITION,
                        y = PERCENT,
                        colour = MASK,
                        group = interaction(MASK, POP, drop = TRUE)
                    )
                ) +
                    geom_line(
                        linewidth = 0.65,
                        alpha = 0.95,
                        linejoin = "round"
                    ) +
                    scale_colour_manual(
                        values = colours,
                        breaks = page_masks,
                        guide = "none"
                    ) +
                    scale_x_continuous(
                        breaks = scales::breaks_pretty(n = 5),
                        expand = expansion(mult = c(0.01, 0.01))
                    ) +
                    scale_y_continuous(
                        limits = c(0, NA),
                        breaks = scales::breaks_pretty(n = 4),
                        expand = expansion(mult = c(0, 0.08))
                    ) +
                    plot_theme +
                    theme(
                        legend.position = "none",
                        panel.spacing = grid::unit(0.8, "lines"),
                        strip.text = element_text(size = 10),
                        strip.text.y = element_text(angle = 0),
                        axis.text.x = element_text(
                            angle = 45,
                            hjust = 1
                        )
                    ) +
                    labs(
                        title = paste(
                            variant_type,
                            if (per_population) {
                                "per-population mask diagnostics:"
                            } else {
                                "mask diagnostics:"
                            },
                            selected_rule
                        ),
                        subtitle = if (per_population) {
                            "Normalised within each mask"
                        } else {
                            "Normalised within each mask"
                        },
                        x = selected_rule,
                        y = "Sites within category (% per bin)"
                    )

                if (per_population) {
                    p + facet_grid(
                        rows = vars(MASK),
                        cols = vars(POP),
                        scales = "fixed",
                        labeller = labeller(
                            MASK = as_labeller(mask_labels),
                            POP = label_value
                        )
                    )
                } else {
                    p + facet_wrap(
                        vars(MASK),
                        ncol = 3,
                        scales = "fixed",
                        labeller = labeller(
                            MASK = as_labeller(mask_labels)
                        )
                    )
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

    # Produce separate files for each observed variant type.
    # REF and ALL are also included if present.
    variant_types <- ordered_types(c(
        as.character(plot_steps$TYPE),
        as.character(record_totals$TYPE)
    ))

    if (!length(variant_types)) {
        stop("No variant types available for plotting.")
    }

    for (variant_type in variant_types) {
        file_type <- tolower(variant_type)

        # Global: one statistic per page, faceted by mask.
        pdf(
            paste0("vcf_masks_", file_type, "_global.pdf"),
            width = 11,
            height = 8
        )

        present_rules <- plot_steps[
            POP == "." & TYPE == variant_type,
            unique(as.character(RULE))
        ]

        rules <- intersect(
            levels(plot_steps$RULE),
            present_rules
        )

        if (length(rules)) {
            for (selected_rule in rules) {
                x <- plot_steps[
                    POP == "." &
                    TYPE == variant_type &
                    as.character(RULE) == selected_rule
                ]

                print(distribution_plot(
                    x,
                    variant_type,
                    selected_rule
                ))

                rm(x)
                invisible(gc())
            }
        } else {
            empty_page(paste(
                "No global mask histogram data for",
                variant_type
            ))
        }

        dev.off()

        # Population: one statistic per page.
        # Masks occupy rows; populations occupy columns.
        pdf(
            paste0("vcf_masks_", file_type, "_perpop.pdf"),
            width = 14,
            height = 10
        )

        present_rules <- plot_steps[
            POP != "." & TYPE == variant_type,
            unique(as.character(RULE))
        ]

        rules <- intersect(
            levels(plot_steps$RULE),
            present_rules
        )

        if (length(rules)) {
            for (selected_rule in rules) {
                x <- plot_steps[
                    POP != "." &
                    TYPE == variant_type &
                    as.character(RULE) == selected_rule
                ]

                print(distribution_plot(
                    x,
                    variant_type,
                    selected_rule,
                    per_population = TRUE
                ))

                rm(x)
                invisible(gc())
            }
        } else {
            empty_page(paste(
                "No per-population mask histogram data for",
                variant_type
            ))
        }

        dev.off()

        # Counts: individual categories plus the masked union.
        pdf(
            paste0("vcf_masks_", file_type, "_counts.pdf"),
            width = 9,
            height = 6
        )

        x <- copy(record_totals[TYPE == variant_type])

        if (nrow(x)) {
            setorder(x, N_RECORDS, MASK)

            x[, LABEL := sub(
                "_FAIL$", "", sub("^MASK_", "", MASK)
            )]

            x[MASK == "AnyMask", LABEL := "Any mask (union)"]

            # Keep the original ID as the axis key.
            x[, MASK := factor(
                MASK,
                levels = as.character(MASK)
            )]

            axis_labels <- setNames(
                x$LABEL,
                as.character(x$MASK)
            )

            p <- ggplot(x, aes(MASK, N_RECORDS)) +
                geom_col(fill = "#619CFF") +
                scale_x_discrete(labels = axis_labels) +
                coord_flip() +
                scale_y_continuous(
                    labels = scales::label_number(big.mark = ","),
                    expand = expansion(mult = c(0, 0.05))
                ) +
                theme_classic() +
                labs(
                    title = paste(variant_type, "mask membership"),
                    x = NULL,
                    y = "Number of sites"
                )

            print(p)
            rm(p)
        } else {
            empty_page(paste(
                "No mask membership data for",
                variant_type
            ))
        }

        dev.off()

        rm(x)
        invisible(gc())
    }
    RSCRIPT
    '''
}