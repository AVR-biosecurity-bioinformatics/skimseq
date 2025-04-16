/*
    Process reads
*/

//// import modules
include { ALIGN_GENOME                    } from '../modules/align_genome'
include { ALIGN_MITO                    } from '../modules/align_mito'
include { CONSENSUS_MITO                         } from '../modules/consensus_mito'
include { FASTP                         } from '../modules/fastp'
include { FASTQC as FASTQC_POSTTRIM     } from '../modules/fastqc'
include { FASTQC as FASTQC_PRETRIM      } from '../modules/fastqc'
include { PROCESS_BAM_MITO                 } from '../modules/process_bam_mito'

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

    ALIGN_MITO (
        FASTP.out.fastq,
        ch_mito_indexed
    )

    // group mito .bam files by sample
    ALIGN_MITO.out.bam
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

    ALIGN_GENOME (
        FASTP.out.fastq,
        ch_genome_indexed
    )



    emit: 
    CONSENSUS_MITO.out.fasta


}