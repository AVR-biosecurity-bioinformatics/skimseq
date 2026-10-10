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
    source(args[[1L]])

    suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
    })
    
    # Reading, validation and aggregation are handled by functions.R.
    mask_data <- read_variant_mask_histograms("histograms")

    # AnyMask is useful for union counts, but redundant in distributions.
    plot_steps <- copy(mask_data$steps[MASK != "AnyMask"])
    record_totals <- copy(mask_data$records)

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

    # Include only categories actually present.
    present_masks <- unique(plot_steps$MASK)

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
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    distribution_plot <- function(x, title, per_population = FALSE) {
        p <- ggplot(
            x,
            aes(
                x = POSITION,
                y = PERCENT,
                colour = MASK,
                group = MASK
            )
        ) +
            geom_step(direction = "hv", linewidth = 0.65) +
            scale_colour_manual(
                values = colours,
                breaks = mask_order,
                labels = unname(mask_labels[mask_order])
            ) +
            scale_x_continuous(
                breaks = scales::breaks_pretty(n = 6)
            ) +
            plot_theme +
            labs(
                title = title,
                subtitle = "Normalised separately within each mask category",
                caption = paste(
                    "Categories overlap. Unmasked does not imply PASS.",
                    "Missing metric values are excluded from denominators."
                ),
                x = NULL,
                y = "Sites within category (% per bin)"
            )

        if (per_population) {
            p + facet_wrap(~POP, scales = "free_y")
        } else {
            p + facet_wrap(~RULE, scales = "free")
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

    # Global distributions: one page per variant type.
    global_df <- plot_steps[POP == "."]

    pdf("vcf_masks_global.pdf", width = 11, height = 8)

    if (nrow(global_df)) {
        for (variant_type in ordered_types(global_df$TYPE)) {
            print(distribution_plot(
                global_df[TYPE == variant_type],
                paste(variant_type, "mask diagnostics")
            ))
        }
    } else {
        empty_page("No global mask histogram data")
    }

    dev.off()

    # Population distributions: one page per metric and variant type.
    per_pop_df <- plot_steps[POP != "."]

    pdf("vcf_masks_perpop.pdf", width = 11, height = 8)

    if (nrow(per_pop_df)) {
        combinations <- unique(per_pop_df[, .(RULE, TYPE)])
        setorder(combinations, TYPE, RULE)

        for (i in seq_len(nrow(combinations))) {
            selected_rule <- as.character(combinations$RULE[[i]])
            variant_type <- combinations$TYPE[[i]]

            x <- per_pop_df[
                as.character(RULE) == selected_rule &
                TYPE == variant_type
            ]

            print(distribution_plot(
                x,
                paste(
                    variant_type,
                    "per-population",
                    selected_rule
                ),
                per_population = TRUE
            ))
        }
    } else {
        empty_page("No per-population mask histogram data")
    }

    dev.off()

    # Counts: individual categories plus the masked union.
    pdf("vcf_masks_counts.pdf", width = 9, height = 6)

    if (nrow(record_totals)) {
        for (variant_type in ordered_types(record_totals$TYPE)) {
            x <- copy(record_totals[TYPE == variant_type])
            setorder(x, N_RECORDS, MASK)

            x[, LABEL := sub(
                "_FAIL$", "", sub("^MASK_", "", MASK)
            )]

            x[MASK == "AnyMask", LABEL := "Any mask (union)"]

            # Use the original ID as the axis key to keep categories distinct.
            x[, MASK := factor(MASK, levels = as.character(MASK))]

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
                    caption = paste(
                        "Individual categories overlap;",
                        "Any mask counts each masked record once."
                    ),
                    x = NULL,
                    y = "Number of sites"
                )

            print(p)
        }
    } else {
        empty_page("No mask membership data")
    }

    dev.off()
    RSCRIPT
    '''
}