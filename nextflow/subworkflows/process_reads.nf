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

workflow PROCESS_READS {

    take:
    ch_reads
    ch_mito_indexed
    ch_genome_indexed

    main: 

    /* 
        Read QC
    */

    FASTQC_PRETRIM (
        ch_reads,
        "pretrim"
    )

    FASTP (
        ch_reads
    )

    // channel for post-trim fastqc
    FASTP.out.fastq
        .map { sample, read1, read2, json -> [ sample, read1, read2 ] }
        .set { ch_fastqc_posttrim_input }

    FASTQC_POSTTRIM (
        ch_fastqc_posttrim_input,
        "posttrim"
    )

    /*
        Mitochondrial variant calling
    */

    // align reads to mitochondrial genome
    MAP_TO_MITO (
        FASTP.out.fastq,
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
        FASTP.out.fastq,
        ch_genome_indexed
    )

    // group nuclear .bam files by sample
    MAP_TO_GENOME.out.bam
        .groupTuple ( by: 0 )
        .set { ch_grouped_genome_bam }

    // process nuclear .bam (merge, sort, index)
    PROCESS_BAM_GENOME (
        ch_grouped_genome_bam
    )

    // extract unmapped reads
    EXTRACT_UNMAPPED (
        PROCESS_BAM_GENOME.out.bam
    )

    // base quality score recalibration (if a list of known variants are provided)


    // generate statistics about the genome .bam files
    BAM_STATS (
        EXTRACT_UNMAPPED.out.bam
    )


    emit: 
    mito_fasta = CONSENSUS_MITO.out.fasta
    bam = BAM_STATS.out.bam
    bam_stats = BAM_STATS.out.stats


}