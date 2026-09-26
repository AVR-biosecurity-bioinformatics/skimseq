/*
    Genotype mitochondrial variants
*/

//// import modules
include { CONSENSUS_MITO                           } from '../modules/consensus_mito/consensus_mito'
include { REALIGN_MITO                             } from '../modules/realign_mito/realign_mito'
include { REALIGN_MITO as REALIGN_MITO_SHIFTED     } from '../modules/realign_mito/realign_mito'
include { PILEUP_MITO                              } from '../modules/pileup_mito/pileup_mito'
include { PILEUP_MITO as PILEUP_MITO_SHIFTED       } from '../modules/pileup_mito/pileup_mito'
include { COUNT_CRAM_PERBASE as COUNT_MITO_PERBASE } from '../modules/count_cram_perbase/count_cram_perbase'
include { COUNT_CRAM_PERBASE as COUNT_MITO_PERBASE_SHIFTED } from '../modules/count_cram_perbase/count_cram_perbase'
include { COMBINE_MOSDEPTH_EVENTS as COMBINE_MITO_EVENTS } from '../modules/combine_mosdepth_events/combine_mosdepth_events'
include { COMBINE_MOSDEPTH_EVENTS as COMBINE_MITO_EVENTS_SHIFTED } from '../modules/combine_mosdepth_events/combine_mosdepth_events'
include { CREATE_INTERVAL_CHUNKS as CREATE_MITO_CHUNKS } from '../modules/create_interval_chunks/create_interval_chunks'
include { CREATE_INTERVAL_CHUNKS as CREATE_MITO_CHUNKS_SHIFTED } from '../modules/create_interval_chunks/create_interval_chunks'
include { GATHER_MITO_PILEUPS as GATHER_MITO_ORIGINAL } from '../modules/gather_mito_pileups/gather_mito_pileups'
include { GATHER_MITO_PILEUPS as GATHER_MITO_SHIFTED  } from '../modules/gather_mito_pileups/gather_mito_pileups'

workflow MITO_GENOTYPING {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_mito_indexed
    ch_shifted_mito_indexed
    ch_mito_bed
    ch_mito_shifted_bed
    ch_numt_bed
    ch_dummy_file
   
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

    REALIGN_MITO.out.mito_cram
        .map { sample, cram, crai ->
            tuple('original', sample, cram, crai)
        }
        .groupTuple(by: 0)
        .set { ch_mito_crams_grouped }

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

    REALIGN_MITO_SHIFTED.out.mito_cram
        .map { sample, cram, crai ->
            tuple('shifted', sample, cram, crai)
        }
        .groupTuple(by: 0)
        .set { ch_shifted_mito_crams_grouped }

    /*
    * Count depth separately on the original and shifted mitochondrial
    * references. Each count process emits tuple(sample, events.starch).
    */
    COUNT_MITO_PERBASE(
        REALIGN_MITO.out.mito_cram,
        ch_mito_indexed,
        ch_dummy_file
    )

    COUNT_MITO_PERBASE_SHIFTED(
        REALIGN_MITO_SHIFTED.out.mito_cram,
        ch_shifted_mito_indexed,
        ch_dummy_file
    )

    // Collect event archives independently for each reference.
    ch_original_mito_events = COUNT_MITO_PERBASE.out.events
        .map { _sample, starch -> starch }
        .toList()
        .filter { archives -> !archives.isEmpty() }

    ch_shifted_mito_events = COUNT_MITO_PERBASE_SHIFTED.out.events
        .map { _sample, starch -> starch }
        .toList()
        .filter { archives -> !archives.isEmpty() }

    // Build a zero-filled cohort RLE track in each reference's
    // own coordinate system.
    COMBINE_MITO_EVENTS(
        ch_mito_indexed,
        ch_mito_bed,
        ch_original_mito_events,
        true
    )

    COMBINE_MITO_EVENTS_SHIFTED(
        ch_shifted_mito_indexed,
        ch_mito_shifted_bed,
        ch_shifted_mito_events,
        true
    )

    // Create workload-balanced intervals independently for each track.
    // The final false allows splitting within the mitochondrial contig.
    CREATE_MITO_CHUNKS(
        ch_mito_bed,
        COMBINE_MITO_EVENTS.out.rle,
        params.mito_bases_per_chunk,
        0,
        false
    )

    CREATE_MITO_CHUNKS_SHIFTED(
        ch_mito_shifted_bed,
        COMBINE_MITO_EVENTS_SHIFTED.out.rle,
        params.mito_bases_per_chunk,
        0,
        false
    )


    /*
     * Generate pileups
     */

    // Pair each original-reference interval with the original-reference cohort.
    CREATE_MITO_CHUNKS.out.interval_bed
        .flatMap { _selector, beds, tbis ->
            def bedList = beds instanceof List ? beds : [beds]
            def tbiList = tbis instanceof List ? tbis : [tbis]

            def tbiByName = tbiList.collectEntries { index ->
                [(index.name): index]
            }

            bedList.collect { bed ->
                def index = tbiByName["${bed.name}.tbi"]

                assert index != null :
                    "Missing index for ${bed.name}"

                tuple(
                    bed.name.replaceFirst(/\.bed\.gz$/, ''),
                    bed,
                    index
                )
            }
        }
        .combine(ch_mito_crams_grouped)
        .map { chunk_id, bed, bed_tbi, label, samples, crams, crais ->
            tuple(label, chunk_id, bed, bed_tbi, samples, crams, crais)
        }
        .set { ch_original_mito_pileup_inputs }

    // Pair each shifted-reference interval with the shifted-reference cohort.
    CREATE_MITO_CHUNKS_SHIFTED.out.interval_bed
        .flatMap { _selector, beds, tbis ->
            def bedList = beds instanceof List ? beds : [beds]
            def tbiList = tbis instanceof List ? tbis : [tbis]

            def tbiByName = tbiList.collectEntries { index ->
                [(index.name): index]
            }

            bedList.collect { bed ->
                def index = tbiByName["${bed.name}.tbi"]

                assert index != null :
                    "Missing index for ${bed.name}"

                tuple(
                    bed.name.replaceFirst(/\.bed\.gz$/, ''),
                    bed,
                    index
                )
            }
        }
        .combine(ch_shifted_mito_crams_grouped)
        .map { chunk_id, bed, bed_tbi, label, samples, crams, crais ->
            tuple(label, chunk_id, bed, bed_tbi, samples, crams, crais)
        }
        .set { ch_shifted_mito_pileup_inputs }

    PILEUP_MITO(
        ch_original_mito_pileup_inputs,
        ch_mito_indexed
    )

    PILEUP_MITO_SHIFTED(
        ch_shifted_mito_pileup_inputs,
        ch_shifted_mito_indexed
    )

     /*
     * Consensus from original and shifted reference pileups
     */

    // Concatenate original-reference counts in chunk order.
    ch_original_counts = PILEUP_MITO.out.counts
        .map { _cohort, chunk_id, _samples_tsv, counts_tsv ->
            tuple(chunk_id, counts_tsv)
        }
        .collectFile(
            name: 'original.all_sites.tsv',
            sort: { entry -> entry[0] }
        ) { entry ->
            entry[1]
        }

    // Concatenate shifted-reference counts independently.
    ch_shifted_counts = PILEUP_MITO_SHIFTED.out.counts
        .map { _cohort, chunk_id, _samples_tsv, counts_tsv ->
            tuple(chunk_id, counts_tsv)
        }
        .collectFile(
            name: 'shifted.all_sites.tsv',
            sort: { entry -> entry[0] }
        ) { entry ->
            entry[1]
        }

    // The corrected PILEUP_MITO writes the same ordered sample
    // manifest for every chunk. Take one copy for consensus.
    ch_original_samples = PILEUP_MITO.out.counts
        .map { _cohort, _chunk_id, samples_tsv, _counts_tsv -> samples_tsv }
        .first()

    ch_original_samples
        .combine(ch_original_counts)
        .combine(ch_shifted_counts)
        .map { samples_tsv, original_counts, shifted_counts ->
            tuple('all', samples_tsv, original_counts, shifted_counts)
        }
        .set { ch_consensus_inputs }

    CONSENSUS_MITO(
        ch_consensus_inputs,
        ch_mito_indexed
    )
    emit: 
    mito_consensus = CONSENSUS_MITO.out.consensus

}