

//// import subworkflows
// include { PROCESS_GENOME                                             } from '../subworkflows/process_genome'
include { FILTER_VARIANTS                                             } from '../subworkflows/filter_variants'
include { GATK_GENOTYPING                                             } from '../subworkflows/gatk_genotyping'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MITO_GENOTYPING                                             } from '../subworkflows/mito_genotyping'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'
include { CREATE_BED_INTERVALS                                              } from '../modules/create_bed_intervals' 
//include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 
//include { CONVERT_INTERVALS                                             } from '../modules/convert_intervals' 

// Import dummny file
ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)

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
    ch_genome_bed = INDEX_GENOME.out.bed.first()

    // Handle empty intervals for interval bed
    if ( params.interval_bed ){
        ch_interval_bed = Channel
            .fromPath (
                 params.interval_bed, 
                 checkIfExists: true
             )
    } else {
        // Set to whole genome
        ch_interval_bed = ch_genome_bed
    } 

    // Handle empty intervals for exclude bed
    if ( params.interval_exclude_bed ){
        ch_exclude_bed = Channel
            .fromPath (
                params.interval_exclude_bed, 
                checkIfExists: true
            )
    } else {
        ch_exclude_bed = ch_dummy_file
    }
    
    // create groups of genomic intervals for parallel genotyping
    CREATE_BED_INTERVALS (
        ch_genome_indexed,
        params.interval_n,
        params.interval_break_n,
        params.interval_break_n_length,
        params.interval_subdivide,
        ch_interval_bed,
        ch_exclude_bed,
        params.interval_exclude_padding,
        params.mito_contig
    )

    // create intervals channel, with one interval_list file per element
    CREATE_BED_INTERVALS.out.interval_list
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
    ch_mito_bed = INDEX_MITO.out.bed.first()
    
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
        ch_mito_bed
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