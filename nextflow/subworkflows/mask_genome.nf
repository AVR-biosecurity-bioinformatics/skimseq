/*
    Mask Genome
*/

//// import modules
include { CREATE_GENOME_MASKS                                       } from '../modules/create_genome_masks' 
include { MERGE_MASKS                                               } from '../modules/merge_masks' 
include { SUMMARISE_MASKS                                           } from '../modules/summarise_masks' 

workflow MASK_GENOME {

    take:
    ch_genome_indexed
    ch_include_bed
    ch_exclude_bed
    ch_mito_bed
    //ch_sample_cram

    main: 

    // Create masks from exclude intervals and existing genome masks
    CREATE_GENOME_MASKS (
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        params.exclude_padding,
        params.use_reference_hardmasks,
        params.use_reference_softmasks
    )

    /*
    Create mask file and summarise
    */

    //Concatenate multiple masks together intp a list
    CREATE_GENOME_MASKS.out.mask_bed
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