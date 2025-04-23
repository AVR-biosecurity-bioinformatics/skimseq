/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                                 } from '../modules/call_variants'
include { CONVERT_INTERVALS                                             } from '../modules/convert_intervals' 
include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 



workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed

    main: 

    // check that params.interval_size is not less than 100k bases
    if ( params.interval_size < 1E5 ){
        error "'params.interval_size' must be >= 1E5 (100,000 bases)"
    }

    /* 
        Read QC
    */

    // create genome intervals for genotyping
    CREATE_INTERVALS (
        ch_genome_indexed,
        params.interval_size
    )

    // turn intervals file into GATK format via .bed
    CONVERT_INTERVALS (
        CREATE_INTERVALS.out.intervals,
        ch_genome_indexed
    )

    // create intervals channel, with one interval_list file per element
    CONVERT_INTERVALS.out.interval_list
        .flatten()
        .set { ch_interval_list }

    // call variants for single samples
    CALL_VARIANTS (
        ch_sample_bam,
        ch_genome_indexed
    )

    // joint genotype GVCFs over genomic intervals

    

    emit: 
    CALL_VARIANTS.out


}