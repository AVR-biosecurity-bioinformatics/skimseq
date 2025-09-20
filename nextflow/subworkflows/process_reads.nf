/*
    Process reads
*/

//// import modules
include { BAM_STATS                             } from '../modules/bam_stats'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome'
include { VALIDATE_FASTQ                        } from '../modules/validate_fastq'
include { SPLIT_FASTQ                           } from '../modules/split_fastq'
include { REPAIR_FASTQ                          } from '../modules/repair_fastq'
include { MERGE_BAM                             } from '../modules/merge_bam'

workflow PROCESS_READS {

    take:
    ch_reads
    ch_genome_indexed

    main: 

    // collect FASTP filtering parameters into a single list
    Channel.of(
        params.rf_quality,
        params.rf_length,
        params.rf_n_bases,
        params.rf_trim_polyg,
        params.rf_cut_right,
        params.rf_cut_window_size,
        params.rf_cut_mean_quality,
        params.rf_lc_filter,
        params.rf_lc_threshold,
        params.rf_correction,
        params.rf_overlap_length,
        params.rf_overlap_diff,
        params.rf_overlap_diff_pc,
        params.rf_custom_flags
    )
    .collect( sort: false )
    .set { ch_fastp_filters }
    

    /* 
        Validate fastq and extract read headers
    */
    VALIDATE_FASTQ (
        ch_reads
    )

    // Convert stdout to a string for status (PASS or FAIL)
    VALIDATE_FASTQ.out.fastq_with_status
        .map { sample, lib, read1, read2, stdout -> [ sample, lib, read1, read2, stdout.trim() ] }
        .branch { sample, lib, read1, read2, status ->
            fail: status == 'FAIL'
            pass: status == 'PASS'
        }
        .set { validation_routes }

    // Print a warning if any samples fail validation and need to be repaired
    validation_routes.fail
    .map { sample, lib, read1, read2, _ -> sample } 
    .unique()
    .collect()
    .map { fails ->
        if (fails && fails.size() > 0)
        log.warn "Repairing malformed FASTQs for ${fails.size()} sample(s): ${fails.join(', ')}"
        true
    }
    .set { _warn_done }  // force evaluation

    // Repair any fastqs that failed validation 
    REPAIR_FASTQ(
        validation_routes.fail.map { sample, lib, read1, read2, _ -> [sample, lib, read1, read2] }
    )

    // Join repaired fastqs back into validated fastqs
    validation_routes.pass.map { sample, lib, read1, read2, _ -> [sample, lib, read1, read2] }
        .mix( REPAIR_FASTQ.out.fastq )
        .set { ch_all_fixed_fastq }

    /* 
        Read splitting
    */

    // split paired fastq files into chunks for parallel processing
    SPLIT_FASTQ (
        ch_all_fixed_fastq,
        params.fastq_chunk_size
    )

    // Create new channel with each fastq chunk
    SPLIT_FASTQ.out.fastq_interval
        .splitCsv ( by: 1, elem: 4, sep: "," )
	    .map { sample, lib, read1, read2, intervals -> [ sample, lib, read1, read2, intervals[0], intervals[1] ] }
	    .set { ch_fastq_split }

    /* 
        Read filtering and alignments
    */

    MAP_TO_GENOME (
        ch_fastq_split,
        ch_fastp_filters,
        ch_genome_indexed
    )
    
    // Merge chunked .bam files by sample (column 0), filter, and index
    MERGE_BAM (
        MAP_TO_GENOME.out.bam.groupTuple ( by: 0 )
    )

    // extract unmapped reads
    // TODO: Make this optional
    EXTRACT_UNMAPPED (
        MERGE_BAM.out.bam
    )

    // TODO: base quality score recalibration (if a list of known variants are provided)

    // generate QC statistics for the merged .bam files
    BAM_STATS (
        MERGE_BAM.out.bam
    )

    // Create reports channel for multiqc
    MAP_TO_GENOME.out.json.map { sample, lib, start, end, json -> [ sample, json ] }
        .mix(BAM_STATS.out.stats, BAM_STATS.out.flagstats, BAM_STATS.out.coverage, MERGE_BAM.out.markdup)
        .set { ch_reports}

    // Create sample renaming table to handle chunks in multiqc report
    MAP_TO_GENOME.out.json
    .map { sample, lib, start, end, json ->
        def unit = "${lib}.${start}-${end}"
        [ unit, sample ]
    }
    .set { renaming_table }
    
    emit: 
    bam = MERGE_BAM.out.bam
    reports = ch_reports
    renaming_table = renaming_table
}

