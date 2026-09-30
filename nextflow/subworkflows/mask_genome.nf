/*
    Mask Genome
*/

//// import modules
include { GENMAP                                                    } from '../modules/genmap/genmap' 
include { LONGDUST                                                  } from '../modules/longdust/longdust'
include { COMBINE_BEDS as COMBINE_MASKS                             } from '../modules/combine_beds/combine_beds' 
include { SUMMARISE_MASKS                                           } from '../modules/summarise_masks/summarise_masks' 

workflow MASK_GENOME {

    take:
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
    - Cohort wide depth variation
    */

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

}