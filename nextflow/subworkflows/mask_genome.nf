nextflow.enable.types = true

/*
 * Mask genome
 */

include { GENMAP                        } from '../modules/genmap/genmap'
include { LONGDUST                      } from '../modules/longdust/longdust'
include { COUNT_CRAM_WINDOWS            } from '../modules/count_cram_windows/count_cram_windows'
include { SPLIT_BED_INTO_WINDOWS        } from '../modules/split_bed_into_windows/split_bed_into_windows'
include { COMBINE_MOSDEPTH_WINDOWS      } from '../modules/combine_mosdepth_windows/combine_mosdepth_windows'
include { COMBINE_BEDS as COMBINE_MASKS } from '../modules/combine_beds/combine_beds'
include { SUMMARISE_MASKS               } from '../modules/summarise_masks/summarise_masks'

// Import record types
include {
    SampleCram
    PreparedGenome
    MaskSummary
} from '../types'

workflow MASK_GENOME {

    take:
    ch_sample_cram: Channel<SampleCram>
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
     * Read-alignment-derived masks
     *
     * - cohort-wide normalised depth
     * - breadth above minimum depth
     */

    windows = SPLIT_BED_INTO_WINDOWS(
        genome
    )

    ch_cram_windows = COUNT_CRAM_WINDOWS(
        ch_sample_cram,
        genome,
        windows
    )

    normalised_windows = ch_cram_windows
        .filter { result ->
            result.normalised_windows != null
        }
        .map { result ->
            result.normalised_windows
        }
        .collect()

    depth_masks = COMBINE_MOSDEPTH_WINDOWS(
        normalised_windows
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
        true,
        0,
        4,
        'distinct'
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
     * Summarise final mask
     */

    mask_summary = SUMMARISE_MASKS(
        genome,
        combined_mask
    )

    emit:
    mask_bed: Value<Path> = combined_mask
    summary: Value<MaskSummary> = mask_summary
}