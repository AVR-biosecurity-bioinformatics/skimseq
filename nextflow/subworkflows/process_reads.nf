/*
    Process reads
*/

//// import modules
include { BAM_STATS                             } from '../modules/bam_stats'
include { DOWNSAMPLE_READS                      } from '../modules/downsample_reads'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { FASTP                                 } from '../modules/fastp'
include { FASTQC as FASTQC_PRETRIM              } from '../modules/fastqc'
include { FASTQTOBAM                            } from '../modules/fastqtobam'
include { SPLIT_FASTQ                           } from '../modules/split_fastq'
//include { FASTQC as FASTQC_POSTTRIM             } from '../modules/fastqc'
//include { PROCESS_BAM_GENOME                    } from '../modules/process_bam_genome'
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
    
    // collect BAM filtering parameters into a single list
    Channel.of(
        params.bam_rmdup
    )
    .collect( sort: false )
    .set { ch_bam_filters }

    /*
        Downsampling and replication of input reads
    */

    Channel.of(params.downsample_reads)
        .map { val -> val.tokenize(',') }
        .flatten()
        .set { ch_downsample_reads }

    if ( !params.downsample_reps.toString().isInteger() ){
        error ("*** PARAMETER ERROR: '--downsample_reps' must be an integer ***")
    }

    Channel.of(1..params.downsample_reps)
        .set { ch_downsample_reps }

    ch_reads
        .combine(ch_downsample_reads)
        .combine(ch_downsample_reps)
        .set { ch_downsample_input }

    if ( ch_downsample_reps ){
        DOWNSAMPLE_READS (
            ch_downsample_input,
            ch_genome_indexed
        )

        ch_reads_ds = DOWNSAMPLE_READS.out.reads

    } else {
        ch_reads_ds = ch_reads
    }

    /* 
        Read QC
    */

    FASTQC_PRETRIM (
        ch_reads_ds,
        "pretrim"
    )

    // split paired fastq into chunks using seqtk
    SPLIT_FASTQ (
        ch_reads_ds,
        params.fastq_chunk_size
    )

    SPLIT_FASTQ.out.fastq_interval
        .splitCsv ( by: 1, elem: 3, sep: "," )
        .map { sample, read1, read2, intervals -> [ sample, read1, read2, intervals[0], intervals[1] ] }
        .set { ch_fastq_split }

    /* 
        Read alignments
    */

    FASTQTOBAM (
        ch_fastq_split,
        ch_fastp_filters,
        ch_genome_indexed,
        ch_bam_filters
    )
    
    // group nuclear .bam files by sample
    FASTQTOBAM.out.bam
        .groupTuple ( by: 0 )
        .set { ch_grouped_genome_bam }

    // extract unmapped reads
    EXTRACT_UNMAPPED (
        ch_grouped_genome_bam
    )

    // TODO: base quality score recalibration (if a list of known variants are provided)

    // generate statistics about the genome .bam files
    BAM_STATS (
        FASTQTOBAM.out.bam
    )

    emit: 
    bam = ch_grouped_genome_bam
    bam_stats = BAM_STATS.out.stats
}