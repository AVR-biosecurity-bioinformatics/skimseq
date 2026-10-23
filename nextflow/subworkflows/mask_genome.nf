nextflow.enable.types = true

/*
 * Mask genome
 */

include { GENMAP                        } from '../modules/genmap/genmap'
include { LONGDUST                      } from '../modules/longdust/longdust'
include { BIGWIG_COHORT_AVERAGE         } from '../modules/bigwig_cohort_average/bigwig_cohort_average'
include { COHORT_DEPTH_MASK             } from '../modules/cohort_depth_mask/cohort_depth_mask'
include { COMBINE_BEDS as COMBINE_MASKS } from '../modules/combine_beds/combine_beds'
include { SUMMARISE_MASKS               } from '../modules/summarise_masks/summarise_masks'

// Import record types
include {
    SampleCram
    CramStats
    PreparedGenome
    MaskSummary
} from '../types'

workflow MASK_GENOME {

    take:
    ch_cram_stats: Channel<CramStats>
    genome: Value<PreparedGenome>
    mito_bed: Value<Path>
    numt_bed: Value<Path>

    main:

    /*
     * Reference-derived masks
     */

    genmap_mask = GENMAP(
        genome
    )

    longdust_mask = LONGDUST(
        genome
    )

    /*
     * Coverage masks
     */

    // Collect all per-sample Riker BigWigs and mean coverage values.
    cohort_cram_stats = ch_cram_stats
        .collect()
        .map { stats ->

            def valid_stats = stats
                .findAll { result ->
                    def keep = result.mean_coverage.toDouble() > 0.0d

                    if (!keep) {
                        log.warn(
                            "Excluding '${result.sample}' from cohort coverage: " +
                            "mean coverage is ${result.mean_coverage}"
                        )
                    }

                    keep
                }
                .toList()

            assert !valid_stats.isEmpty() :
                "No samples have positive mean coverage"

            record(
                samples: valid_stats
                    .collect { result -> result.sample }
                    .toList(),

                bigwigs: valid_stats
                    .collect { result -> result.wgs_per_base }
                    .toList(),

                mean_coverages: valid_stats
                    .collect { result -> result.mean_coverage }
                    .toList()
            )
        }

    // Scale each BigWig by 1 / mean coverage and calculate the
    // cohort-wide mean normalised coverage in fixed-size windows.
    depth_masks = BIGWIG_COHORT_AVERAGE(
        cohort_cram_stats
    )

    // Test figwig version
    cohort_bigwigs = ch_cram_stats
        .collect()
        .map { stats ->
            stats.collect { result ->
                result.wgs_per_base
            }
        }

    depth_masks_test = COHORT_DEPTH_MASK(
        cohort_bigwigs
    )

    /*
     * Combine masks
     */

    reference_mask = genome
        .map { result ->
            result.reference_masks
        }

    low_depth_mask = depth_masks
        .map { result ->
            result.low_mask
        }

    high_depth_mask = depth_masks
        .map { result ->
            result.high_mask
        }

    mask_paths = channel
        .of(0)
        .combine(reference_mask)
        .combine(genmap_mask)
        .combine(longdust_mask)
        .combine(low_depth_mask)
        .combine(high_depth_mask)
        .combine(numt_bed)
        .combine(mito_bed)
        .map { masks ->
            [
                masks[1],
                masks[2],
                masks[3],
                masks[4],
                masks[5],
                masks[6],
                masks[7]
            ]
        }

    combined_mask_channel = COMBINE_MASKS(
        mask_paths,
        genome,
        false, // Preserve individual labelled intervals.
        0,     // Ignored when merging is disabled.
        4,     // Ignored.
        'distinct' // Ignored.
    )

    combined_mask = combined_mask_channel
        .collect()
        .map { masks ->

            if (masks.size() != 1) {
                error(
                    "Expected exactly one combined mask BED, " +
                    "found ${masks.size()}"
                )
            }

            masks.toList()[0]
        }

    /*
     * Summarise final masks
     */

    mask_summary = SUMMARISE_MASKS(
        genome,
        combined_mask
    )

    emit:
    mask_bed: Value<Path> = combined_mask
    summary: Value<MaskSummary> = mask_summary
}