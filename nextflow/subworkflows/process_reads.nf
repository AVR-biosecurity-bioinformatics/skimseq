/*
    Process reads
*/

//// import modules
include { BAM_STATS                             } from '../modules/bam_stats'
include { MAP_TO_GENOME                         } from '../modules/map_to_genome'
include { MAP_TO_MITO                           } from '../modules/map_to_mito'
include { CONSENSUS_MITO                        } from '../modules/consensus_mito'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { FASTP                                 } from '../modules/fastp'
include { FASTQC as FASTQC_POSTTRIM             } from '../modules/fastqc'
include { FASTQC as FASTQC_PRETRIM              } from '../modules/fastqc'
include { PROCESS_BAM_GENOME                    } from '../modules/process_bam_genome'
include { PROCESS_BAM_MITO                      } from '../modules/process_bam_mito'
include { SPLIT_FASTQ                           } from '../modules/split_fastq'

workflow PROCESS_READS {

    take:
    ch_reads
    ch_mito_indexed
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
        Read QC
    */

    FASTQC_PRETRIM (
        ch_reads,
        "pretrim"
    )

    FASTP (
        ch_reads,
        ch_fastp_filters
    )

    // channel for post-trim fastqc
    FASTP.out.fastq
        .map { sample, read1, read2, json -> [ sample, read1, read2 ] }
        .set { ch_fastqc_posttrim_input }

    FASTQC_POSTTRIM (
        ch_fastqc_posttrim_input,
        "posttrim"
    )

    // split paired fastq into chunks using seqtk
    SPLIT_FASTQ (
        FASTP.out.fastq,
        params.fastq_chunk_size
    )

    // combine matching chunks from paired files
    SPLIT_FASTQ.out.fastq
        .transpose()
        .multiMap { sample, file1, file2, json ->
            first: [ sample, json, file1 ]
            second:  [ sample, json, file2 ]
        }
        .set { ch_split_multi }

    ch_split_multi.first
        .map { sample, json, reads_file ->
            def filename_list = reads_file.getFileName().toString().split("\\.")
            [ sample, json, filename_list[1], reads_file ]
        }
        .set { ch_split_read1 }
    
    ch_split_multi.second
        .map { sample, json, reads_file ->
            def filename_list = reads_file.getFileName().toString().split("\\.")
            [ sample, json, filename_list[1], reads_file ]
        }
        .set { ch_split_read2 }

    ch_split_read1
        .join ( 
            ch_split_read2, 
            by: [0,1,2],
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { sample, json, chunk, file1, file2 ->
            [ sample, file1, file2, json ]
        }
        .set { ch_fastq_split }


    /*
        Mitochondrial variant calling
    */

    // align reads to mitochondrial genome
    MAP_TO_MITO (
        ch_fastq_split,
        ch_mito_indexed
    )

    // group mito .bam files by sample
    MAP_TO_MITO.out.bam
        .groupTuple ( by: 0 )
        .set { ch_grouped_mito_bam }

    // process mito bam (merge, sort, index)
    PROCESS_BAM_MITO (
        ch_grouped_mito_bam
    )

    // call consensus fasta file from mito bam
    CONSENSUS_MITO (
        PROCESS_BAM_MITO.out.bam,
        ch_mito_indexed
    )

    /* 
        Nuclear variant calling
    */

    // align reads to nuclear genome 
    MAP_TO_GENOME (
        ch_fastq_split,
        ch_genome_indexed
    )

    // group nuclear .bam files by sample
    MAP_TO_GENOME.out.bam
        .groupTuple ( by: 0 )
        .set { ch_grouped_genome_bam }

    // process nuclear .bam (merge, sort, index)
    PROCESS_BAM_GENOME (
        ch_grouped_genome_bam,
        ch_bam_filters
    )

    // extract unmapped reads
    EXTRACT_UNMAPPED (
        PROCESS_BAM_GENOME.out.bam
    )

    // base quality score recalibration (if a list of known variants are provided)


    // generate statistics about the genome .bam files
    BAM_STATS (
        PROCESS_BAM_GENOME.out.bam
    )


    emit: 
    mito_fasta = CONSENSUS_MITO.out.fasta
    bam = PROCESS_BAM_GENOME.out.bam
    bam_stats = BAM_STATS.out.stats


}