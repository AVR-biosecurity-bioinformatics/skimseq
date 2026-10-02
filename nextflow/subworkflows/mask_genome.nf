/*
    Mask Genome
*/

//// import modules
include { GENMAP                            } from '../modules/genmap/genmap' 
include { LONGDUST                          } from '../modules/longdust/longdust'
include { COUNT_CRAM_WINDOWS                } from '../modules/count_cram_windows/count_cram_windows'
include { SPLIT_BED_INTO_WINDOWS            } from '../modules/split_bed_into_windows/split_bed_into_windows'
include { COMBINE_MOSDEPTH_WINDOWS          } from '../modules/combine_mosdepth_windows/combine_mosdepth_windows'
include { COMBINE_BEDS as COMBINE_MASKS     } from '../modules/combine_beds/combine_beds' 
include { SUMMARISE_MASKS                   } from '../modules/summarise_masks/summarise_masks' 

workflow MASK_GENOME {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_calling_bed
    ch_reference_masks
    ch_mito_bed
    ch_numt_bed

    main: 

    // Create mapabillity mask with GENMAP
    GENMAP (
       ch_genome_indexed,
       params.genmap_kmer_length,
       params.genmap_error_tol,
       params.genmap_thresh
    )

    // Create Repeat/LCR mask with longdust
    LONGDUST (
       ch_genome_indexed,
       params.longdust_kmer_length,
       params.longdust_window_size,
       params.longdust_thresh
    )

    /*
    Read alignment based masks
    - Cohort wide normalised depth
    - breadth above minimum depth
    */
    
    // TODO: Change 100 to a prameter
    SPLIT_BED_INTO_WINDOWS(
        ch_calling_bed, 
        ch_genome_indexed,
        1000
     )

    // Count per-base read depths in all crams, used for masking
    // TODO: Could later output perbase from this and use for sample coverage calculates
    COUNT_CRAM_WINDOWS(
        ch_sample_cram,
        ch_genome_indexed,
        SPLIT_BED_INTO_WINDOWS.out.windows.first()
    )

    ch_normalised_windows = COUNT_CRAM_WINDOWS.out.normalised_windows
        .map { _sample, bed, _tbi -> bed }
        .collect()

    COMBINE_MOSDEPTH_WINDOWS(ch_normalised_windows)

    /*
    Variant calling based masks
    - Cohort wide allele balance
    -
    */

    /*
    Create final file and summarise
    */

    //Concatenate multiple masks together intp a list
    ch_reference_masks
      .concat(GENMAP.out.mask_bed)
      .concat(LONGDUST.out.mask_bed)
      .concat(COMBINE_MOSDEPTH_WINDOWS.out.low_mask)
      .concat(COMBINE_MOSDEPTH_WINDOWS.out.high_mask)
      .concat(ch_numt_bed)
      .concat(ch_mito_bed)
      .collect()
      .set{ ch_mask_bed }

    // Merge all masks
    COMBINE_MASKS (
        ch_mask_bed,
        ch_genome_indexed,
        true,
        0,
        4,
        "distinct"
    )

    // Summarise masks
    SUMMARISE_MASKS (
        ch_genome_indexed,
        ch_calling_bed,
        COMBINE_MASKS.out.bed
    )

    emit: 
    mask_bed = COMBINE_MASKS.out.bed
    mask_summary = SUMMARISE_MASKS.out.mask_summary
    mask_summary_bed = SUMMARISE_MASKS.out.mask_summary_bed
    mask_pass_bed = SUMMARISE_MASKS.out.mask_pass_bed
    perbase = channel.empty()
}