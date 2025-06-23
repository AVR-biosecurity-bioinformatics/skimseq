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

    // process mito bam (merge, sort, index)
    PROCESS_BAM_MITO (
        ch_sample_bam,
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