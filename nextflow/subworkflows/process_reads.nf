/*
    Process reads
*/

//// import modules
include { BAM_STATS                             } from '../modules/bam_stats'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome'
include { SPLIT_FASTQ                           } from '../modules/split_fastq'
include { MERGE_BAM                             } from '../modules/merge_bam'

//include { FASTQC as FASTQC_POSTTRIM             } from '../modules/fastqc'
//include { MAP_TO_GENOME                         } from '../modules/map_to_genome'

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


    /* 
        Read splitting
    */

    // split paired fastq into chunks using seqtk
    SPLIT_FASTQ (
        ch_reads,
        params.fastq_chunk_size
    )

    SPLIT_FASTQ.out.fastq_interval
        .splitCsv ( by: 1, elem: 3, sep: "," )
	    .map { sample, read1, read2, intervals -> [ sample, read1, read2, intervals[0], intervals[1] ] }
	    .set { ch_fastq_split }

    /* 
        Read filtering and alignments
    */

    MAP_TO_GENOME (
        ch_fastq_split,
        ch_fastp_filters,
        ch_genome_indexed
    )
    
    // group chunked .bam files by sample
    MAP_TO_GENOME.out.bam
        .groupTuple ( by: 0 )
        .set { ch_grouped_genome_bam }

    // Merge chunked .bam files by sample, filter, and index
    MERGE_BAM (
        ch_grouped_genome_bam
    )

    // extract unmapped reads
    EXTRACT_UNMAPPED (
        MERGE_BAM.out.bam
    )

    // TODO: base quality score recalibration (if a list of known variants are provided)

    // generate statistics about the .bam files
    BAM_STATS (
        MERGE_BAM.out.bam
    )

    // Create reports channel for multiqc
    MAP_TO_GENOME.out.json
        .mix(BAM_STATS.out.stats, BAM_STATS.out.flagstats, BAM_STATS.out.coverage, MERGE_BAM.out.markdup)
        .set { ch_reports}

    emit: 
    bam = MERGE_BAM.out.bam
    reports = ch_reports
}