nextflow.enable.types = true

/*
    Genotype mitochondrial variants
*/

//// import modules
include { NUMT_MASK                                              } from '../modules/numt_mask/numt_mask'
include { REALIGN_MITO                                           } from '../modules/realign_mito/realign_mito'
include { REALIGN_MITO as REALIGN_MITO_SHIFTED                   } from '../modules/realign_mito/realign_mito'
include { PILEUP_MITO                                            } from '../modules/pileup_mito/pileup_mito'
include { PILEUP_MITO as PILEUP_MITO_SHIFTED                     } from '../modules/pileup_mito/pileup_mito'
include { CREATE_INTERVAL_CHUNKS as CREATE_MITO_CHUNKS           } from '../modules/create_interval_chunks/create_interval_chunks'
include { CREATE_INTERVAL_CHUNKS as CREATE_MITO_CHUNKS_SHIFTED   } from '../modules/create_interval_chunks/create_interval_chunks'
include { CONSENSUS_MITO                                         } from '../modules/consensus_mito/consensus_mito'

// Import record types
include {
    SampleCram
    PreparedGenome
    PreparedMito
    MitoConsensus
} from '../types'

workflow MITO_GENOTYPING {

    take:
    sample_cram: Channel<SampleCram>
    genome: Value<PreparedGenome>
    mito: Value<PreparedMito>
    cohort_size: Value<Integer>

    main: 

    /*
    * Calculate the number of mitochondrial chunks.
    */

    mito_included_bases = mito
        .map { prepared ->
            prepared.included_bases.text.trim().toLong()
        }

    mito_n_chunks = mito_included_bases
        .combine(cohort_size)
        .map { included_bases_value, cohort_size_value ->

            def total_bases =
                cohort_size_value * included_bases_value

            def bases_per_chunk =
                params.mito_target_sample_bases

            Math.max(
                1,
                (
                    (total_bases + bases_per_chunk - 1) /
                    bases_per_chunk
                ) as int
            )
        }


    /*
    * Identify putative NUMT regions in the nuclear genome.
    */

    numt_bed = NUMT_MASK(
        genome,
        mito
    )

    /*
    * Realign reads to the original and shifted mitochondrial references.
    *
    * Both calls recruit reads using the original mitochondrial interval,
    * but use different alignment targets.
    */

    original_mito_crams = REALIGN_MITO(
        sample_cram,
        genome,
        mito,
        numt_bed,
        false
    )

    shifted_mito_crams = REALIGN_MITO_SHIFTED(
        sample_cram,
        genome,
        mito,
        numt_bed,
        true
    )


    /*
    * Create workload-balanced mitochondrial chunks.
    */

    // Beds need to be extracted from the mito record for compatibility with existing chunk code
    original_chunks = CREATE_MITO_CHUNKS(
        mito.map { prepared_mito -> prepared_mito.mito_bed},
        mito_n_chunks,
        true // split large intervals (required for mito)
    )

    shifted_chunks = CREATE_MITO_CHUNKS_SHIFTED(
        mito.map { prepared_mito -> prepared_mito.shifted_mito_bed },
        mito_n_chunks,
        true // split large intervals (required for mito)
    )

    /*
    * Expand the original-reference chunk bundle into one record per chunk.
    */

    original_intervals = original_chunks
        .flatMap { chunks ->
            chunks.beds.collect { bed ->

                def index = chunks.tbis.find { candidate ->
                    candidate.name == "${bed.name}.tbi"
                }

                if (index == null) {
                    throw new IllegalStateException(
                        "Missing tabix index for ${bed.name}"
                    )
                }

                record(
                    reference: 'original',
                    interval_id: bed.name.replaceFirst(/\.bed\.gz$/, ''),
                    interval_bed: bed,
                    interval_index: index
                )
            }
        }


    /*
    * Expand the shifted-reference chunk bundle into one record per chunk.
    */

    shifted_intervals = shifted_chunks
        .flatMap { chunks ->
            chunks.beds.collect { bed ->

                def index = chunks.tbis.find { candidate ->
                    candidate.name == "${bed.name}.tbi"
                }

                if (index == null) {
                    throw new IllegalStateException(
                        "Missing tabix index for ${bed.name}"
                    )
                }

                record(
                    reference: 'shifted',
                    interval_id: bed.name.replaceFirst(/\.bed\.gz$/, ''),
                    interval_bed: bed,
                    interval_index: index
                )
            }
        }


    /*
    * Append the corresponding cohort to each interval record.
    */

    original_cohort = original_mito_crams.collect()
    shifted_cohort = shifted_mito_crams.collect()

    original_mito_pileup_inputs = original_intervals
        .combine(samples: original_cohort)

    shifted_mito_pileup_inputs = shifted_intervals
        .combine(samples: shifted_cohort)

    /*
    * Generate pileups.
    */

    original_pileups = PILEUP_MITO(
        original_mito_pileup_inputs,
        mito,
        false // dont use shifted input
    )

    shifted_pileups = PILEUP_MITO_SHIFTED(
        shifted_mito_pileup_inputs,
        mito,
        true // use shifted input
    )

    /*
    * Collect pileup outputs for consensus.
    */

    original_pileups_collected = original_pileups
        .collect()
        .map { pileups ->
            pileups.toSorted { pileup -> pileup.interval_id }
        }

    shifted_pileups_collected = shifted_pileups
        .collect()
        .map { pileups ->
            pileups.toSorted { pileup -> pileup.interval_id }
        }

    consensus_input = original_pileups_collected
        .map { pileups ->
            def first_pileup = pileups.find { _pileup -> true }

            if (first_pileup == null) {
                throw new IllegalStateException(
                    "No original mitochondrial pileups were produced"
                )
            }

            record(
                samples_tsv: first_pileup.samples_tsv,
                original_counts: pileups.collect { pileup ->
                    pileup.counts_tsv
                }
            )
        }
        .combine(
            shifted_counts: shifted_pileups_collected.map { pileups ->
                pileups.collect { pileup ->
                    pileup.counts_tsv
                }
            }
        )

    mito_consensus = CONSENSUS_MITO(
        consensus_input,
        mito
    )

    emit:
    mito_consensus: Value<MitoConsensus> = mito_consensus
    numt_bed: Value<Path> = numt_bed
}