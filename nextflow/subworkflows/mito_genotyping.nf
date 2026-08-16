/*
    Genotype mitochondrial variants
*/

//// import modules
include { CONSENSUS_MITO                        } from '../modules/consensus_mito/consensus_mito'
include { REALIGN_MITO                          } from '../modules/realign_mito/realign_mito'
include { PILEUP_MITO                           } from '../modules/pileup_mito/pileup_mito'


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

    // Pileup mitochondrial reads (whole cohort)
    REALIGN_MITO.out.mito_bams
        .map { sample, bam, bai, shifted_bam, shifted_bai -> tuple('all', sample, bam, bai, shifted_bam, shifted_bai ) }
        .groupTuple(by: 0)
        .set{ ch_mito_bams_grouped }

    PILEUP_MITO(
        ch_mito_bams_grouped,
        ch_mito_indexed,
        ch_shifted_mito_indexed
    )

    // Call consensus from pileup
    CONSENSUS_MITO(
        PILEUP_MITO.out.counts,
        ch_mito_indexed
    )

    // Align consensus mito reads

    emit: 
    mito_consensus = CONSENSUS_MITO.out.consensus

}