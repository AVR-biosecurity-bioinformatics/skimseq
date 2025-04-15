/*
    Process reads
*/

//// import modules
include { FASTQC as FASTQC_PRETRIM } from '../modules/fastqc'

workflow PROCESS_READS {

    take:
    ch_reads

    main: 


    FASTQC_PRETRIM (
        ch_reads,
        "pretrim"
    )


    emit: 
    FASTQC_PRETRIM.out


}