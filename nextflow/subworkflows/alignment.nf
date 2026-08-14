/*
    Process reads
*/

//// import modules
include { VALIDATE_CRAM                         } from '../modules/validate_cram/validate_cram'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome/map_to_genome'
include { SPLIT_FASTQ                           } from '../modules/split_fastq/split_fastq'
include { MERGE_CRAM                            } from '../modules/merge_cram/merge_cram'
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
        Find and validate any pre-existing crams, these will be skipped
        To pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    */

    // Use existing crams if they are present and the option is set
    //if( params.use_existing_cram ) {
    //    ch_sample_names
    //        .map { sample ->
    //            def cram = file("${params.cram_store}/${sample}.cram")
    //            def crai = file("${cram}.crai")
    //            tuple(sample, cram, crai)
    //        }
    //        .filter { sample, cram, crai -> cram.exists() && crai.exists() }
    //        .set { ch_existing_cram }
    //
    //    // Validate cram files by default
    //    if( !params.skip_cram_validation ) {
    //        VALIDATE_CRAM (
    //            ch_reads.join(ch_existing_cram, by: 0),
    //            ch_genome_indexed
    //        )
    //
    //        // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
    //        VALIDATE_CRAM.out.status
    //            .map { sample, stdout -> tuple(sample, stdout.trim()) }
    //            .join( ch_existing_cram, by: 0 )
    //            .map { sample, status, cram, crai -> tuple(sample, cram, crai, status) }
    //            .branch {  sample, cram, crai, status ->
    //                fail: status == 'FAIL'
    //                pass: status == 'PASS'
    //                invalid: true
    //            }
    //            .set { cram_validation_routes }
    //
    //        // Fail loudly if there is an invalid status
    //        cram_validation_routes.invalid
    //            .map { sample, cram, crai, status ->
    //                throw new IllegalStateException(
    //                    "Unexpected CRAM validation status for ${sample}: '${status}'"
    //                )
    //            }
    //            .set { _invalid_cram_status }
    //
    //        // Channel with just passing crams
    //        cram_validation_routes.pass
    //            .map { sample, cram, crai, status -> tuple(sample, cram, crai) } 
    //            .set { ch_validated_cram }
    //            
    //        // Print warning if any cram files exist but fail validation
    //        cram_validation_routes.fail
    //            .map { sample, cram, crai, status -> sample }
    //            .unique()
    //            .collect()
    //            .subscribe { fails ->
    //                if (fails) {
    //                    log.warn(
    //                        "CRAM validation failed for ${fails.size()} sample(s): " +
    //                        fails.join(', ')
    //                    )
    //                }
    //            }
    //
    //    } else {
    //      // Skip validation, assume all existing crams are good
    //      ch_validated_cram = ch_existing_cram 
    //    }
    //
    //    // Sample ids that already have a good CRAM
    //    ch_validated_cram
    //        .map { sample, cram, crai -> sample }
    //        .toList()
    //        .map { ids -> ids as Set } 
    //        .set { ch_cram_done }
    //} else{
    //    ch_cram_done = Channel.value([] as Set)
    //    ch_validated_cram = channel.empty()
    //}

    // Temporary - remap all crams
    ch_cram_done = Channel.value([] as Set)
    ch_validated_cram = channel.empty()
    
    // Filter the reads to only those samples who dont already have a validated cram - only these will be mapped
    ch_reads
        .combine(ch_cram_done)  
        .filter { sample, lib, source, read1, read2, local_reads, doneSet -> !(doneSet as Set).contains(sample) }
        .map {  sample, lib, source, read1, read2, local_reads, doneSet -> tuple(sample, lib, source, read1, read2, local_reads) }
        .set { ch_reads_to_map }

    // Group all FASTQ pairs and libraries belonging to each sample.
    ch_reads_to_map
        .groupTuple(by: 0)
        .map {
            sample,
            libs,
            sources,
            input1s,
            input2s,
            local_reads_groups ->

            def unique_sources = sources.unique()

            if (unique_sources.size() != 1) {
                error(
                    "Sample '${sample}' contains multiple input source types: " +
                    "${unique_sources.join(', ')}. Mixed source types are not " +
                    "currently supported within one mapping task."
                )
            }

            def source = unique_sources.first()

            def grouped_local_reads = source == 'local'
                ? local_reads_groups.flatten()
                : []

            tuple(
                sample,
                libs,
                source,
                input1s,
                input2s,
                grouped_local_reads
            )
        }
        .set { ch_reads_grouped_by_sample }

    /*
    * Prepare mapping intervals.
    * Each channel emission represents one complete library, potentially containing multiple FASTQ pairs.
    * n_intervals is the number of libraries that must finish before merging
    * Sorted by remote, then largest to smallest local
    * This ensures longer running jobs start first, smaller jobs backfill once resources available
    */
    ch_reads_grouped_by_lib
        .groupTuple(by: 0)
        .flatMap {
            sample,
            libs,
            sources,
            input1_groups,
            input2_groups,
            local_reads_groups ->

            int n_intervals = libs.size()

            /*
            * Total size of all local FASTQs belonging to this sample.
            * Remote libraries contribute zero.
            */
            long sample_local_size = (0..<n_intervals).sum { i ->
                sources[i] == 'local'
                    ? local_reads_groups[i].sum { read ->
                        read.size()
                    } as long
                    : 0L
            } as long

            (0..<n_intervals).collect { i ->
                boolean is_local = sources[i] == 'local'

                /*
                * Temporary sorting tuple:
                *
                *   0: priority: remote=0, local=1
                *   1: total local FASTQ size for the sample
                *   2: sample
                *   3: number of libraries for the sample
                *   4: library
                *   5: source
                *   6: list of R1 inputs for the library
                *   7: list of R2 inputs for the library
                *   8: flattened local FASTQs for the library
                */
                tuple(
                    is_local ? 1 : 0,
                    sample_local_size,
                    sample,
                    n_intervals,
                    libs[i],
                    sources[i],
                    input1_groups[i],
                    input2_groups[i],
                    local_reads_groups[i]
                )
            }
        }
        .toSortedList { a, b ->
            (a[0] <=> b[0]) ?:  // Remote libraries first
            (b[1] <=> a[1]) ?:  // Largest local samples first
            (a[2] <=> b[2]) ?:  // Keep each sample together
            (a[4] <=> b[4])     // Deterministic library order
        }
        .flatMap { it }
        .map { priority, sample_local_size, sample, n_intervals, lib, source, input1s, input2s, local_reads ->
            tuple(
                sample,
                n_intervals,
                lib,
                source,
                input1s,
                input2s,
                local_reads
            )
        }
        .set { ch_reads_to_map_intervals }
    /* 
        Read mapping
    */

    // Align reads to genome
    MAP_TO_GENOME (
        ch_reads_to_map_intervals,
        ch_genome_indexed
    )

    // Print warning if any files had different numbers of forward and reverse reads
    MAP_TO_GENOME.out.fastq_warnings
        .map { sample, lib, warning_file ->
            tuple(sample, lib, warning_file.text.trim())
        }
        .subscribe { sample, lib, warning ->
            log.warn("MAP_TO_GENOME ${sample}:${lib}: ${warning}")
        }
    
    // Collect mapping CRAMs as soon as all FASTQ pairs for a sample complete
    MAP_TO_GENOME.out.cram
        .map { sample, n_intervals, lib, cram, crai ->
            tuple(groupKey(sample, n_intervals), cram, crai)
        }
        .groupTuple()
        .map { key, crams, crais ->
            tuple(key.getGroupTarget(), crams, crais)
        }
        .set { ch_cram_to_merge }

    // Merge chunked .cram files by sample
    MERGE_CRAM(
        ch_cram_to_merge,
        ch_genome_indexed
    )

    // combine validated existing CRAMs with newly created CRAMs
    ch_validated_cram
      .mix( MERGE_CRAM.out.cram )
      .distinct { it[0] }      // dedupe by sample if needed
      .set{ ch_sample_cram }

    // Helper process to stage intermediate CRAMs 
    STAGE_CRAM(
        ch_sample_cram
    )

    // Count per-base depths in cram, used for masking and creating interval chunks
    COUNT_CRAM_PERBASE (
        STAGE_CRAM.out.cram,
        ch_genome_indexed,
        ch_exclude_bed
    )

    emit: 
    cram = STAGE_CRAM.out.cram
    perbase = COUNT_CRAM_PERBASE.out.perbase
    counts = COUNT_CRAM_PERBASE.out.counts

}

