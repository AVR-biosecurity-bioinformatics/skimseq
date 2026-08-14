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
    // Sort by remote then largest to smallest local files.
    // Slower alignment tasks run first, quick jobs backfill
    ch_reads_to_map
        .groupTuple(by: 0)
        .map {
            sample, libs, sources, input1s, input2s, local_reads_groups ->
            
            // Check only one input source is defined per process
            // TODO? Could allow multiple input sources in future
            def unique_sources = sources.unique(false)
            if (unique_sources.size() != 1) {
                error(
                    "Sample '${sample}' contains multiple input source " +
                    "types: ${unique_sources.join(', ')}. Mixed local, URL " +
                    "and accession inputs are not currently supported " +
                    "within one MAP_TO_GENOME task."
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
                    "R1=${input1s.size()}, " +
                    "R2=${input2s.size()}, " +
                    "local read groups=${local_reads_groups.size()}."
                )
            }

            def source = unique_sources.first()

            /*
             * Each local samplesheet row contributes [R1, R2].
             * Flatten to: row1_R1, row1_R2, row2_R1, row2_R2, ...
             * MAP_TO_GENOME can reconstruct these using collate(2).
             */
            def grouped_local_reads = source == 'local'
                ? local_reads_groups.flatten()
                : []

            long sample_local_size = source == 'local'
                ? grouped_local_reads.sum { read ->
                    read.size()
                } as long
                : 0L

            tuple(
                source == 'local' ? 1 : 0, // priority: remote=0, local=1
                sample_local_size,
                sample,
                libs,
                source,
                input1s,
                input2s,
                grouped_local_reads
            )
        }
        .toSortedList { a, b ->
            (a[0] <=> b[0]) ?:  // Remote samples first
            (b[1] <=> a[1]) ?:  // Largest local samples first
            (a[2] <=> b[2])     // Deterministic sample order
        }
        .flatMap { sorted_samples ->
            sorted_samples
        }
        .map {
            priority, sample_local_size, sample, libs, source, input1s, input2s, local_reads ->
            tuple( sample, libs, source, input1s, input2s, local_reads )
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

    emit: 
    cram = STAGE_CRAM.out.cram
    perbase = COUNT_CRAM_PERBASE.out.perbase
    counts = COUNT_CRAM_PERBASE.out.counts

}

