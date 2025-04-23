/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                             } from '../modules/call_variants'
include { CONVERT_INTERVALS                                              } from '../modules/convert_intervals' 
include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 



workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed

    main: 

    /* 
        Read QC
    */

    // create genome intervals for genotyping
    CREATE_INTERVALS (
        ch_genome_indexed
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

    

    emit: 
    CALL_VARIANTS.out


}