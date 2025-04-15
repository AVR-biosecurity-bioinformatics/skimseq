/*
    Process reads
*/

//// import modules
include { FASTP                     } from '../modules/fastp'
include { FASTQC as FASTQC_PRETRIM } from '../modules/fastqc'

workflow PROCESS_READS {

    take:
    ch_reads

    main: 


    FASTQC_PRETRIM (
        ch_reads,
        "pretrim"
    )

    FASTP (
        ch_reads
    )



    emit: 
    FASTQC_PRETRIM.out


}