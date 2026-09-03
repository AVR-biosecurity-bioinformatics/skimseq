/*
    Genotype mitochondrial variants
*/

//// import modules
include { CONSENSUS_MITO                        } from '../modules/consensus_mito/consensus_mito'
include { REALIGN_MITO                          } from '../modules/realign_mito/realign_mito'
include { REALIGN_MITO as REALIGN_MITO_SHIFTED  } from '../modules/realign_mito/realign_mito'
include { PILEUP_MITO                           } from '../modules/pileup_mito/pileup_mito'
include { PILEUP_MITO as PILEUP_MITO_SHIFTED    } from '../modules/pileup_mito/pileup_mito'


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
     * Realign to original mito reference
     */

    REALIGN_MITO(
        ch_sample_cram,
        ch_genome_indexed,
        ch_mito_indexed,
        ch_mito_bed,
        ch_numt_bed
    )

    REALIGN_MITO.out.mito_bam
        .map { sample, bam, bai ->
            tuple('all', sample, bam, bai)
        }
        .groupTuple(by: 0)
        .set { ch_mito_bams_grouped }

    /*
     * Realign to shifted mito reference
     */
    REALIGN_MITO_SHIFTED(
        ch_sample_cram,
        ch_genome_indexed,
        ch_shifted_mito_indexed,
        ch_mito_bed,
        ch_numt_bed
    )

    REALIGN_MITO_SHIFTED.out.shifted_mito_bam
        .map { sample, bam, bai ->
            tuple('all', sample, bam, bai)
        }
        .groupTuple(by: 0)
        .set { ch_shifted_mito_bams_grouped }

    /*
     * Generate pileups
     */
    PILEUP_MITO(
        ch_mito_bams_grouped,
        ch_mito_indexed
    )

    PILEUP_MITO_SHIFTED(
        ch_shifted_mito_bams_grouped,
        ch_shifted_mito_indexed
    )

    /*
     * Consensus from original and shifted reference pileups
     */
    PILEUP_MITO.out.counts
        .join(PILEUP_MITO_SHIFTED.out.counts, by: 0 )
        .map {cohort, samples_tsv, original_counts, _shifted_samples_tsv, shifted_counts ->
            tuple(cohort, samples_tsv, original_counts, shifted_counts )
        }
        .set { ch_consensus_inputs }

    CONSENSUS_MITO(
        ch_consensus_inputs,
        ch_mito_indexed
    )

    emit: 
    mito_consensus = CONSENSUS_MITO.out.consensus

}