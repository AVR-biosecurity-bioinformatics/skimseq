

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { FILTER_VARIANTS                                             } from '../subworkflows/filter_variants'
include { GATK_GENOTYPING                                             } from '../subworkflows/gatk_genotyping'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MITO_GENOTYPING                                             } from '../subworkflows/mito_genotyping'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'
include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 
include { CONVERT_INTERVALS                                             } from '../modules/convert_intervals' 

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

    //if ( params.mito_genome ){
    //    ch_mito = Channel
    //        .fromPath (
    //            params.mito_genome, 
    //            checkIfExists: true
    //        )
    //} else {
    //    ch_mito = Channel.empty()
    //} 

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
    Process nuclear genome and create intervals
    */

    INDEX_GENOME (
        ch_genome
    )

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()
    
    // TODO: Replace interval creation with bedtools
    // TODO: If bedfile is provided, use that
    
    // create genome intervals for genotyping
    CREATE_INTERVALS (
        ch_genome_indexed,
        params.interval_size
    )

    // split intervals file into chunks of 50 lines for conversion
    CREATE_INTERVALS.out.intervals
        .splitText ( by: 50, file: true )
        .set { ch_intervals }

    // turn intervals into GATK format via .bed
    CONVERT_INTERVALS (
        ch_intervals,
        ch_genome_indexed
    )

    // create intervals channel, with one interval_list file per element
    CONVERT_INTERVALS.out.interval_list
        .flatten()
        // get hash from interval_list name as element to identify intervals
        .map { interval_list ->
            def interval_hash = interval_list.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_list ] }
        .set { ch_interval_list }
        
    /*
    Process mitochondrial genome and create intervals
    */
        
    // TODO: Extract mitochondrial genome contig from reference contig
    INDEX_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = INDEX_MITO.out.fasta_indexed.first()

    // PROCESS_GENOME (
    //     "dummy"
    // )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_reads,
        ch_genome_indexed
    )

    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        PROCESS_READS.out.bam,
        ch_mito_indexed,
        INDEX_MITO.out.bed
    )
    
    /*
    Call variants per sample, then combine and joint-genotype across genomic intervals
    */

    GATK_GENOTYPING (
        PROCESS_READS.out.bam,
        ch_genome_indexed,
        ch_interval_list
    )

    /*
    Filter SNP and INDEL variants
    */

    FILTER_VARIANTS (
        GATK_GENOTYPING.out.vcf,
        ch_genome_indexed
    )

}