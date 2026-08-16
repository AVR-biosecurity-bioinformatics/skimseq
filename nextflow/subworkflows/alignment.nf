/*
    Process reads
*/

//// import modules
include { VALIDATE_CRAM                         } from '../modules/validate_cram/validate_cram'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome/map_to_genome'
include { STAGE_CRAM                            } from '../modules/stage_cram/stage_cram'
include { COUNT_CRAM_PERBASE                    } from '../modules/count_cram_perbase/count_cram_perbase'

workflow ALIGNMENT {

    take:
    ch_sample_names
    ch_reads
    ch_genome_indexed
    ch_exclude_bed

    main: 
    
    /*
    * Group all input rows belonging to each sample.
    */
    ch_reads
        .groupTuple(by: 0)
        .map { sample, libs, sources, input1s, input2s, local_reads_groups ->

            def unique_sources = sources.unique(false)
            if (unique_sources.size() != 1) {
                error(
                    "Sample '${sample}' contains multiple input source types: " +
                    "${unique_sources.join(', ')}. Mixed local, URL, and " +
                    "accession inputs are not currently supported within one " +
                    "MAP_TO_GENOME task."
                )
            }

            if (
                libs.size() != sources.size() ||
                libs.size() != input1s.size() ||
                libs.size() != input2s.size() ||
                libs.size() != local_reads_groups.size()
            ) {
                error(
                    "Input metadata is inconsistent for sample '${sample}': " +
                    "libs=${libs.size()}, " +
                    "sources=${sources.size()}, " +
                    "input1=${input1s.size()}, " +
                    "input2=${input2s.size()}, " +
                    "local read groups=${local_reads_groups.size()}."
                )
            }

            def source = unique_sources.first()

            if (source == 'local') {
                local_reads_groups.eachWithIndex { pair, i ->
                    if (!(pair instanceof Collection) || pair.size() != 2) {
                        error(
                            "Invalid local FASTQ pair for sample '${sample}', " +
                            "row ${i + 1}: ${pair}. Expected [R1, R2]."
                        )
                    }
                }
            }

            def local_r1s = source == 'local'
                ? local_reads_groups.collect { pair -> pair[0] }
                : []

            def local_r2s = source == 'local'
                ? local_reads_groups.collect { pair -> pair[1] }
                : []

            // Return tuple
            tuple(sample, libs, source, input1s, input2s, local_r1s, local_r2s )
        }
        .set { ch_reads_grouped }

    /* 
        Find and validate any pre-existing crams, these will be skipped
        To pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    */
    
    if (params.use_existing_cram) {

        ch_sample_names
            .map { sample ->
                def cram = file("${params.cram_store}/${sample}.cram")
                def crai = file("${cram}.crai")
                tuple(sample, cram, crai)
            }
            .filter { sample, cram, crai -> cram.exists() && crai.exists() }
            .set { ch_existing_cram }

        if (!params.skip_cram_validation) {

            //Join grouped reads to candidate CRAMs.
            ch_reads_grouped
                .join(ch_existing_cram, by: 0)
                .set { ch_cram_validation_input }

            VALIDATE_CRAM(
                ch_cram_validation_input,
                ch_genome_indexed
            )

            // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
            VALIDATE_CRAM.out.status
                .map { sample, stdout -> tuple(sample, stdout.trim()) }
                .join(ch_existing_cram, by: 0)
                .map { sample, status, cram, crai -> tuple(sample, cram, crai, status) }
                .branch { sample, cram, crai, status ->
                    fail: status == 'FAIL'
                    pass: status == 'PASS'
                    invalid: true
                }
                .set { cram_validation_routes }

            // Fail loudly if there is an invalid status
            cram_validation_routes.invalid
                .map { sample, cram, crai, status ->
                    throw new IllegalStateException(
                        "Unexpected CRAM validation status for ${sample}: '${status}'"
                    )
                }
                .set { _invalid_cram_status }

            //Compatible existing CRAMs.
            cram_validation_routes.pass
                .map { sample, cram, crai, status -> tuple(sample, cram, crai) }
                .set { ch_validated_cram }

            // Report incompatible CRAMs. Thes get remapped
            cram_validation_routes.fail
                .map { sample, cram, crai, status -> sample }
                .unique()
                .collect()
                .subscribe { failed_samples ->
                    if (failed_samples) {
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

        //Set of sample names that do not need mapping.
        ch_validated_cram
            .map { sample, cram, crai -> sample }
            .toList()
            .map { ids -> ids as Set } 
            .set { ch_cram_done }

    } else {
        ch_cram_done = channel.value([] as Set)
        ch_validated_cram = channel.empty()
    }


    // Filter the reads to only those samples who dont already have a validated cram - only these will be mapped
    ch_reads_grouped
        .combine(ch_cram_done)
        .filter { sample, libs, source, input1s, input2s, local_r1s, local_r2s, done_set -> !(done_set as Set).contains(sample)}
        .map {sample, libs, source, input1s, input2s, local_r1s, local_r2s, done_set -> tuple(sample, libs, source, input1s, input2s, local_r1s, local_r2s) }
        .set { ch_reads_to_map }

    /*
    * Run remote samples first, followed by the largest local samples.
    */
    ch_reads_to_map
        .map { sample, libs, source, input1s, input2s, local_r1s, local_r2s ->
            long sample_local_size = source == 'local'
                ? (local_r1s + local_r2s).sum { read ->
                    read.size()
                } as long
                : 0L

            tuple(
                source == 'local' ? 1 : 0,
                sample_local_size,
                sample,
                libs,
                source,
                input1s,
                input2s,
                local_r1s,
                local_r2s
            )
        }
        .toSortedList { a, b ->
            (a[0] <=> b[0]) ?:   // Remote first
            (b[1] <=> a[1]) ?:   // Largest local samples first
            (a[2] <=> b[2])      // Deterministic sample order
        }
        .flatMap { sorted_samples ->
            sorted_samples
        }
        .map { priority, sample_local_size, sample, libs, source, input1s, input2s, local_r1s, local_r2s ->
            tuple( sample, libs, source, input1s, input2s,local_r1s, local_r2s )
        }
        .set { ch_reads_grouped_by_sample }

    /* 
        Read mapping
    */

    // Align reads to genome, input is all libraries and reads per sample
    // Output is sample-level cram, no merging required
    MAP_TO_GENOME (
        ch_reads_grouped_by_sample,
        ch_genome_indexed
    )

    // Print warning if any files had different numbers of forward and reverse reads
    MAP_TO_GENOME.out.fastq_warnings
        .map { sample, warning_file ->
            tuple(
                sample,
                warning_file.text.trim()
            )
        }
        .subscribe { sample, warning ->
            log.warn(
                "MAP_TO_GENOME ${sample}: ${warning}"
            )
        }
    
    // Combine pre-validated crams with newly mapped crams
    ch_validated_cram
        .mix(MAP_TO_GENOME.out.cram)
        .distinct { sample, cram, crai -> sample }
        .set { ch_sample_cram }

    // Helper process to stage intermediate CRAMs 
    STAGE_CRAM(
        ch_sample_cram
    )

    // Count per-base depths in all crams, used for masking and creating interval chunks
    COUNT_CRAM_PERBASE (
        STAGE_CRAM.out.cram,
        ch_genome_indexed,
        ch_exclude_bed
    )

    // Only newly generated CRAMs should be published.
    MAP_TO_GENOME.out.cram
        .set { ch_new_cram }

    emit: 
    cram = STAGE_CRAM.out.cram
    new_cram = ch_new_cram
    perbase = COUNT_CRAM_PERBASE.out.perbase
    counts = COUNT_CRAM_PERBASE.out.counts

}

