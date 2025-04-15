

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'


//// import modules



workflow SKIMSEQ {

    /*
    Input channel parsing
    */    

    ch_reads = Channel
        .fromFilePairs(
            "./test/*.R{1,2}.fastq",
            checkIfExists: true, 
            flat: true
        )

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