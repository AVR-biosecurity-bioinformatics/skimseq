/*
    Mask Genome
*/

//// import modules
include { EXTRACT_GENOME_MASKS                                      } from '../modules/extract_genome_masks' 
include { GENMAP                                                    } from '../modules/genmap' 
include { LONGDUST                                                  } from '../modules/longdust'
include { MERGE_MASKS                                               } from '../modules/merge_masks' 
include { SUMMARISE_MASKS                                           } from '../modules/summarise_masks' 

workflow MASK_GENOME {

    take:
    ch_genome_indexed
    ch_include_bed
    ch_exclude_bed
    ch_mito_bed
    ch_read_counts

    main: 

    // Create masks from exclude intervals and existing genome masks
    EXTRACT_GENOME_MASKS (
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        params.exclude_padding,
        params.use_reference_hardmasks,
        params.use_reference_softmasks
    )

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

    // TODO: Create read depth masks from CRAM counts

    /*
    Create mask file and summarise
    */

    //Concatenate multiple masks together intp a list
    EXTRACT_GENOME_MASKS.out.mask_bed
      .concat(GENMAP.out.mask_bed)
      .concat(LONGDUST.out.mask_bed)
      .concat(ch_mito_bed)
      .collect()
      .set{ ch_mask_bed }

    // Merge all masks
    MERGE_MASKS (
        ch_mask_bed
    )
    
    // Summarise masks
    SUMMARISE_MASKS (
        ch_genome_indexed,
        ch_include_bed,
        MERGE_MASKS.out.merged_masks
    )
    
    emit: 
    mask_bed = MERGE_MASKS.out.merged_masks

}