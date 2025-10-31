/*
    Process reads
*/

//// import modules
include { MAP_TO_GENOME                         } from '../modules/map_to_genome'
include { SPLIT_FASTQ                           } from '../modules/split_fastq'
include { MERGE_CRAM                            } from '../modules/merge_cram'

workflow PROCESS_READS {

    take:
    ch_reads_to_map
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
        Read splitting
    */

    // split paired fastq files into chunks for parallel processing
    SPLIT_FASTQ (
        ch_reads_to_map.map { sample, lib, fcid, lane, platform, read1, read2 -> [ sample, lib, read1, read2 ] },
        params.fastq_chunk_size
    )

    // Create new channel with each fastq chunk
    SPLIT_FASTQ.out.fastq_interval
        .splitCsv ( by: 1, elem: 2, sep: "," )
        .map { sample, lib, intervals -> [ sample, lib, intervals[0], intervals[1] ] }
        .join(ch_reads_to_map, by:[0,1] )
        .map { sample, lib, int1, int2, fcid, lane, platform, read1, read2 -> [ sample, lib, fcid, lane, platform, read1, read2, int1, int2 ] }
        .set { ch_fastq_split }

    /* 
        Read filtering and alignments
    */

    MAP_TO_GENOME (
        ch_fastq_split,
        ch_fastp_filters,
        ch_genome_indexed
    )
    
    // Merge chunked .cram files by sample (column 0), filter, and index
    MERGE_CRAM (
        MAP_TO_GENOME.out.cram.groupTuple ( by: 0 ),
        ch_genome_indexed
    )

    // TODO: base quality score recalibration (if a list of known variants are provided)

    // Create sample renaming table to handle chunks in multiqc report
    MAP_TO_GENOME.out.json
    .map { sample, lib, start, end, json ->
        def unit = "${lib}.${start}-${end}"
        [ unit, sample ]
    }
    .set { renaming_table }
    
    emit: 
    cram = MERGE_CRAM.out.cram
    reports = MAP_TO_GENOME.out.json.map { sample, lib, start, end, json -> [ sample, json ] }
    renaming_table = renaming_table
}

