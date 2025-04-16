/*
    Process reads
*/

//// import modules
include { ALIGN_MITO                     } from '../modules/align_mito'
include { FASTP                     } from '../modules/fastp'
include { FASTQC as FASTQC_POSTTRIM } from '../modules/fastqc'
include { FASTQC as FASTQC_PRETRIM } from '../modules/fastqc'

workflow PROCESS_READS {

    take:
    ch_reads
    ch_mito_indexed
    ch_genome

    main: 


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

    ALIGN_MITO (
        FASTP.out.fastq,
        ch_mito_indexed
    )

    // group mito .bam files by sample


    // sort mito .bam, merging if necessary


    // index mito .bam


    // ALIGN_GENOME (

    // )



    emit: 
    FASTQC_PRETRIM.out


}