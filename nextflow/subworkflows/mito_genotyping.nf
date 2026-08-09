/*
    Genotype mitochondrial variants
*/

//// import modules
include { CONSENSUS_MITO                        } from '../modules/consensus_mito/consensus_mito'
include { REALIGN_MITO                          } from '../modules/realign_mito/realign_mito'


workflow MITO_GENOTYPING {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_mito_indexed
    ch_shifted_mito_indexed
    ch_mito_bed
    ch_numt_bed
   

    main: 

    /*
        Mitochondrial variant calling
    */

    // Extract mitochondrial & reads from genomic cram and realign
    REALIGN_MITO (
        ch_sample_cram,
        ch_genome_indexed,
        ch_mito_indexed,
        ch_shifted_mito_indexed,
        ch_mito_bed,
        ch_numt_bed
    )

    // call consensus fasta file from mito bam
    CONSENSUS_MITO (
        ch_sample_cram,
        ch_genome_indexed,
        ch_mito_indexed,
        ch_mito_bed,
        ch_numt_bed,
        params.mito_min_vaf,
        params.mito_min_depth
    )

    // Align consensus mito reads

    emit: 
    mito_fasta = CONSENSUS_MITO.out.fasta

}