

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'


//// import modules



workflow SKIMSEQ {

    /*
    Input channel parsing
    */    

    ch_reads = Channel.fromPath("./test/*.fastq")

    /*
    Process genome 
    */

    // PROCESS_GENOME (
    //     "dummy"
    // )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_reads
    )


}