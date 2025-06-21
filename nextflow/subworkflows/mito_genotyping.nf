/*
    Genotype mitochondrial variants
*/

//// import modules
include { MAP_TO_MITO                           } from '../modules/map_to_mito'
include { CONSENSUS_MITO                        } from '../modules/consensus_mito'
include { PROCESS_BAM_MITO                      } from '../modules/process_bam_mito'


workflow MITO_GENOTYPING {

    take:
    ch_sample_bam
    ch_mito_indexed
    ch_mito_bed

    main: 

    /*
        Mitochondrial variant calling
    */
    // group mito .bam files by sample
    ch_sample_bam
            .flatMap { sample, bam_files, bam_index_files ->
                // For each sample, return individual BAM files with their corresponding BAM index
                bam_files.collect { bam ->
                    def bam_index = "${bam}.bai"  // Assuming the index file has the same name as the BAM file with a .bai extension
                    return tuple(sample, bam, bam_index)  // Return a tuple with sample, BAM file, and BAM index
                }
            }
     .set { ch_grouped_mito_bam }

    // process mito bam (merge, sort, index)
    PROCESS_BAM_MITO (
        ch_grouped_mito_bam,
        ch_mito_bed
    )

    // call consensus fasta file from mito bam
    CONSENSUS_MITO (
        PROCESS_BAM_MITO.out.bam,
        ch_mito_indexed
    )

    emit: 
    mito_fasta = CONSENSUS_MITO.out.fasta

}