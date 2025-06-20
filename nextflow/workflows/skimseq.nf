

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { FILTER_VARIANTS                                             } from '../subworkflows/filter_variants'
include { GATK_GENOTYPING                                             } from '../subworkflows/gatk_genotyping'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'


//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'


workflow SKIMSEQ {

    /*
    Input channel parsing
    */    

    if ( params.samplesheet ){
        ch_samplesheet = Channel
            .fromPath (
                params.samplesheet,
                checkIfExists: true
            )
    } else {
        println "\n*** ERROR: 'params.samplesheet' must be given ***\n"
    }
    
    ch_samplesheet 
        .splitCsv ( by: 1, skip: 1 )
        .map { row -> [ 
            row[0],                                 // sample
            file( row[1], checkIfExists: true ),    // read1 
            file( row[2], checkIfExists: true )     // read2
            ] }
        .set { ch_reads }

    if ( params.mito_genome ){
        ch_mito = Channel
            .fromPath (
                params.mito_genome, 
                checkIfExists: true
            )
    } else {
        ch_mito = Channel.empty()
    } 

    if ( params.ref_genome ){
        ch_genome = Channel
            .fromPath (
                params.ref_genome, 
                checkIfExists: true
            )
    } else {
        ch_genome = Channel.empty()
    } 


    /*
    Process nuclear and mitochondrial genome files
    */

    INDEX_MITO (
        ch_mito
    )

    INDEX_GENOME (
        ch_genome
    )

    // PROCESS_GENOME (
    //     "dummy"
    // )

    ch_mito_indexed = INDEX_MITO.out.fasta_indexed.first()

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_reads,
        ch_mito_indexed,
        ch_genome_indexed
    )

    /*
    Call variants per sample, then combine and joint-genotype across genomic intervals
    */

    GATK_GENOTYPING (
        PROCESS_READS.out.bam,
        ch_genome_indexed
    )

    /*
    Filter SNP and INDEL variants
    */

    //FILTER_VARIANTS (
    //    GATK_GENOTYPING.out.vcf,
    //    ch_genome_indexed
    //)

}