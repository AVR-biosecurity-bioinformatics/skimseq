nextflow.enable.types = true

/*
    Align reads to genome
*/

//// import modules
include { VALIDATE_CRAM                         } from '../modules/validate_cram/validate_cram'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome/map_to_genome'
include { STAGE_CRAM                            } from '../modules/stage_cram/stage_cram'
include { CRAM_STATS_RIKER                      } from '../modules/cram_stats_riker/cram_stats_riker'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped/extract_unmapped'


// Import record types
include {
    SampleReads
    PreparedGenome
    SampleCram
    CramStats
    UnmappedReads
} from '../types'


workflow ALIGNMENT {

    take:
    ch_sample_names: Channel<String>
    ch_reads_grouped: Channel<SampleReads>
    ch_genome_indexed: Value<PreparedGenome>

    main: 

    /* 
        Find and validate any pre-existing crams, these will be skipped from alignment
        To pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    */
    
    if (params.use_existing_cram) {

        ch_existing_cram = ch_sample_names
            .map { sample ->
                def cram = file("${params.cram_store}/${sample}.cram")
                def crai = file("${cram}.crai")

                record(
                    sample: sample,
                    cram: cram,
                    crai: crai
                )
            }
            .filter { candidate ->
                candidate.cram.exists() &&
                    candidate.crai.exists()
            }

        if (!params.skip_cram_validation) {

            //Join grouped reads to candidate CRAMs.
            ch_cram_validation_input = ch_reads_grouped.join(
                ch_existing_cram,
                by: 'sample'
            )

            ch_cram_validation_results = VALIDATE_CRAM(
                ch_cram_validation_input,
                ch_genome_indexed
            )
                .map { result ->
                    def status = result.status.trim()

                    if (status != 'PASS' && status != 'FAIL') {
                        throw new IllegalStateException(
                            "Unexpected CRAM validation status for " +
                                "${result.sample}: '${status}'"
                        )
                    }

                    record(
                        sample: result.sample,
                        cram: result.cram,
                        crai: result.crai,
                        status: status
                    )
                }

            //Compatible existing CRAMs.
            // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
            ch_validated_cram = ch_cram_validation_results
                .filter { result ->
                    result.status == 'PASS'
                }
                .map { result ->
                    record(
                        sample: result.sample,
                        cram: result.cram,
                        crai: result.crai
                    )
                }

            // Report incompatible CRAMs. These get remapped.
            ch_cram_validation_results
                .filter { result ->
                    result.status == 'FAIL'
                }
                .map { result ->
                    result.sample
                }
                .unique()
                .collect()
                .subscribe { failed_samples ->
                    if (!failed_samples.isEmpty()) {
                        log.warn(
                            "CRAM validation failed for " +
                                "${failed_samples.size()} sample(s): " +
                                failed_samples.join(', ')
                        )
                    }
                }

        } else {
            // Validation explicitly disabled: assume all discovered CRAMs are compatible.
            ch_validated_cram = ch_existing_cram
        }

    } else {
        ch_validated_cram = channel.empty()
    }

    ch_cram_done = ch_validated_cram
        .map { sample_cram ->
            sample_cram.sample
        }
        .collect()
    /* 
        Pre-process and align reads to genome
        This is handled in single mapping step
    */

    // Dont map samples which already have a validated CRAM.
    ch_reads_to_map = ch_reads_grouped
        .combine(done_samples: ch_cram_done)
        .filter { reads ->
            !reads.done_samples.contains(reads.sample)
        }
        .map { reads ->
            record(
                sample: reads.sample,
                libs: reads.libs,
                source: reads.source,
                input1s: reads.input1s,
                input2s: reads.input2s,
                local_r1s: reads.local_r1s,
                local_r2s: reads.local_r2s
            )
        }

    // Map samples which dont have validated cram
    mapping_results = MAP_TO_GENOME(
        ch_reads_to_map,
        ch_genome_indexed
    )

    ch_new_cram = mapping_results

    // Combine pre-validated crams with newly generated crams
    ch_sample_cram = ch_validated_cram
        .mix(ch_new_cram)
        .unique { sample_cram ->
            sample_cram.sample
        }

    // Helper process to stage intermediate CRAMs - Handles resume issues
    ch_staged_cram = STAGE_CRAM( ch_sample_cram )

    // generate QC statistics for the merged .cram files
    ch_cram_stats = CRAM_STATS_RIKER(
        ch_staged_cram,
        ch_genome_indexed
    )

    // Optional: extract unmapped reads
    ch_unmapped_reads = channel.empty()

    if (params.output_unmapped_reads) {
        ch_unmapped_reads = EXTRACT_UNMAPPED(
            ch_staged_cram,
            ch_genome_indexed
        )
    }

    emit:
    cram: Channel<SampleCram> = ch_staged_cram
    new_cram: Channel<SampleCram> = ch_new_cram
    cram_stats: Channel<CramStats> = ch_cram_stats
    unmapped_reads: Channel<UnmappedReads> = ch_unmapped_reads
}

