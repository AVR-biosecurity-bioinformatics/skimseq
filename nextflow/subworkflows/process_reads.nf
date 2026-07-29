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

workflow PROCESS_READS {

    take:
    ch_sample_names
    ch_reads
    ch_rg_to_validate
    ch_genome_indexed    

    main: 
    
    /* 
        Find and validate any pre-existing crams, these will be skipped
        To pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    */

    // Use existing crams if they are present and the option is set
    if( params.use_existing_cram ) {
        ch_sample_names
            .map { sample ->
                def cram = file("${params.cram_store}/${sample}.cram")
                def crai = file("${cram}.crai")
                tuple(sample, cram, crai)
            }
            .filter { sample, cram, crai -> cram.exists() && crai.exists() }
            .set { ch_existing_cram }

        // Validate cram files by default
        if( !params.skip_cram_validation ) {
            VALIDATE_CRAM (
                ch_rg_to_validate.join(ch_existing_cram, by: 0),
                ch_genome_indexed
            )

            // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
            VALIDATE_CRAM.out.status
                .map { sample, stdout -> [ sample, stdout.trim() ] }
                .join( ch_existing_cram, by: 0 )
                .map { sample, status, cram, crai -> [ sample, cram, crai, status ] }
                .branch {  sample, cram, crai, status ->
                    fail: status == 'FAIL'
                    pass: status == 'PASS'
                }
                .set { cram_validation_routes }

            // Channel with just passing crams
            cram_validation_routes.pass
                .map { sample, cram, crai, status -> [ sample, cram, crai ] } 
                .set { ch_validated_cram }
                
            // Print warning if any cram files exist but fail validation
            cram_validation_routes.fail
                .map {  sample, cram, crai, status -> sample } 
                .unique()
                .collect()
                .map { fails ->
                    if (fails && fails.size() > 0)
                    log.warn "CRAM file failed validation for ${fails.size()} samples(s): ${fails.join(', ')}"
                    true
                }
                .set { _warn_cram_done }  // force evaluation

        } else {
          // Skip validation, assume all existing crams are good
          ch_validated_cram = ch_existing_cram 
        }

        // Sample ids that already have a good CRAM
        ch_validated_cram
            .map { sample, cram, crai -> sample }
            .toList()
            .map { ids -> ids as Set } 
            .set { ch_cram_done }
    } else{
        ch_cram_done = Channel.value([] as Set)
        ch_validated_cram = channel.empty()
    }

    // Filter the reads to only those samples who dont already have a validated cram - only these will be mapped
    ch_reads
        .combine(ch_cram_done)  
        .filter { sample, lib, fcid, lane, platform, read1, read2, doneSet -> !(doneSet as Set).contains(sample) }
        .map { sample, lib, fcid, lane, platform, read1, read2, doneSet -> tuple(sample, lib, fcid, lane, platform, read1, read2) }
        .set { ch_reads_to_map }

    /* 
        Read splitting
    */

    // Split paired fastq files into even chunks for parallel processing
    SPLIT_FASTQ (
        ch_reads_to_map.map { sample, lib, fcid, lane, platform, read1, read2 -> [ sample, lib, read1, read2 ] },
        params.fastq_chunk_size
    )

    // Parse per-library chunk counts, then sum to per-sample total
    SPLIT_FASTQ.out.nchunks
        .map { sample, lib, nchunks_file ->
            tuple(sample, nchunks_file.text.trim().toInteger())
        }
        .groupTuple(by: 0)
        .map { sample, counts ->
            tuple(sample, counts.sum())
        }
        .set { ch_sample_nchunks }

    // Create new channel with each fastq chunk
    SPLIT_FASTQ.out.fastq_interval
        .splitCsv(by: 1, elem: 2, sep: ",")
        .map { sample, lib, intervals ->
            tuple(sample, lib, intervals[0], intervals[1])
        }
        .combine(ch_reads_to_map, by: [0,1])
        .map { sample, lib, int1, int2, fcid, lane, platform, read1, read2 ->
            tuple(sample, lib, fcid, lane, platform, read1, read2, int1, int2)
        }
        .combine(ch_sample_nchunks, by: 0)
        .map { sample, lib, fcid, lane, platform, read1, read2, int1, int2, n_chunks ->
            tuple(sample, n_chunks, lib, fcid, lane, platform, read1, read2, int1, int2)
        }
        .set { ch_fastq_split }

    /* 
        Read mapping
    */

    // Align reads to genome
    MAP_TO_GENOME (
        ch_fastq_split,
        ch_genome_indexed
    )
    
    // Grouping by sample, nchunks allows early per-sample merge rather than waiting for all MAP_TO_GENOME to finish
    MAP_TO_GENOME.out.cram
        // Add expected group size so each sample emits once all interval CRAms are complete
        .map { sample, n_chunks, cram, crai ->
            tuple(groupKey(sample, n_chunks), cram, crai)
        }
        // Group interval CRAMs by sample, emitting early when n_intervals have arrived
        .groupTuple()
        // Emit sample with grouped crams and indexes for concatenation
        .map { key, crams, crais ->
            tuple(key.getGroupTarget(), crams, crais)
        }
        .set { ch_cram_to_merge }

    // Merge chunked .cram files by sample
    MERGE_CRAM(
        ch_cram_to_merge,
        ch_genome_indexed
    )

    // TODO: base quality score recalibration (if a list of known variants are provided)

    // combine validated existing CRAMs with newly created CRAMs
    ch_validated_cram
      .mix( MERGE_CRAM.out.cram )
      .distinct { it[0] }      // dedupe by sample if needed
      .set{ ch_sample_cram }

    // Helper process to publish CRAMs to output directory. 
    // NOTE: This process (using deep caching) is necessary to avoid violating cache of later steps when inputs switch to existing cram on resume
    STAGE_CRAM(
        ch_sample_cram
    )

    // Count per-base depths in cram, used for masking and creating interval chunks
    COUNT_CRAM_PERBASE (
        STAGE_CRAM.out.cram,
        ch_genome_indexed,
        params.rmdup,
        params.minbq,
        params.minmq
    )

    emit: 
    cram = STAGE_CRAM.out.cram
    perbase = COUNT_CRAM_PERBASE.out.perbase

}

