/*
    Genotype mitochondrial variants
*/

//// import modules
include { MAP_TO_MITO                           } from '../modules/map_to_mito'
include { CONSENSUS_MITO                        } from '../modules/consensus_mito'
include { PROCESS_CRAM_MITO                     } from '../modules/process_cram_mito'


workflow MITO_GENOTYPING {

    take:
    ch_sample_cram
    ch_mito_indexed
    ch_mito_bed
    ch_genome_indexed

    main: 

    /*
        Mitochondrial variant calling
    */

    // Extract mitochondrial reads from genomic cram
    PROCESS_CRAM_MITO (
        ch_sample_cram,
        ch_mito_bed,
        ch_genome_indexed
    )

    // call consensus fasta file from mito bam
    CONSENSUS_MITO (
        PROCESS_CRAM_MITO.out.bam,
        ch_mito_indexed
    )

    emit: 
    mito_fasta = CONSENSUS_MITO.out.fasta

}