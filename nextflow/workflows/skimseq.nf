

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'


//// import modules
// include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'


workflow SKIMSEQ {

    /*
    Input channel parsing
    */    

    ch_reads = Channel
        .fromFilePairs(
            "./test/*_R{1,2}_*.fastq",
            checkIfExists: true, 
            flat: true
        )

    // ch_reads.view()

    if ( params.mito_genome ){
        ch_mito = Channel
            .fromPath(
                params.mito_genome, 
                checkIfExists: true
            )
    } else {
        ch_mito = Channel.empty()
    } 

    if ( params.ref_genome ){
        ch_genome = Channel
            .fromPath(
                params.ref_genome, 
                checkIfExists: true
            )
    } else {
        ch_genome = Channel.empty()
    } 


    /*
    Process genome 
    */

    INDEX_MITO (
        ch_mito
    )

    // INDEX_GENOME (
    //     ch_genome
    // )

    // PROCESS_GENOME (
    //     "dummy"
    // )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_reads,
        INDEX_MITO.out.fasta_indexed,
        ch_genome
    )


}