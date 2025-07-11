/*
    Mask Genome
*/

//// import modules
include { CREATE_GENOME_MASKS                                       } from '../modules/create_genome_masks' 
include { CREATE_GENOME_BINS                                        } from '../modules/create_genome_bins'
include { COUNT_READS_BINS                                          } from '../modules/count_reads_bins'
include { CREATE_MASKS_BINS                                         } from '../modules/create_masks_bins'
include { MERGE_MASKS                                               } from '../modules/merge_masks' 
include { SUMMARISE_MASKS                                           } from '../modules/summarise_masks' 

workflow MASK_GENOME {

    take:
    ch_genome_indexed
    ch_include_bed
    ch_exclude_bed
    ch_mito_bed
    ch_sample_bam

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
    Calculate coverages
    */
    
    // First divide the genome into bins for calculating coverage
    CREATE_GENOME_BINS (
        ch_genome_indexed,
        ch_include_bed,
        params.bin_size
    )

    ch_binned_bed = CREATE_GENOME_BINS.out.binned_bed.first()
    ch_annot_bins = CREATE_GENOME_BINS.out.annotated_bins.first()

    // Count reads in each group of binned intervals
    COUNT_READS_BINS (
          ch_sample_bam,
          ch_binned_bed,
          ch_genome_indexed
    ) 

    // collect counts.tsvs into a single element
    COUNT_READS_BINS.out.counts
        .collect()
        .set { ch_bin_counts }
        
    // Run filter counts module
    CREATE_MASKS_BINS (
          ch_bin_counts,
          ch_binned_bed,
          ch_annot_bins,
          ch_genome_indexed,
          params.bin_gc_lower,
          params.bin_gc_upper,
          params.bin_min_reads,
          params.bin_lower_read_perc,
          params.bin_upper_read_perc,
          params.bin_filter_perc_samples
    )   

    /*
    Create mask file and summarise
    */

    //Concatenate multiple masks together intp a list
    CREATE_GENOME_MASKS.out.mask_bed
      .concat(ch_mito_bed, CREATE_MASKS_BINS.out.bin_masked)
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