/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                                 } from '../modules/call_variants'
include { COMBINE_GVCFS                                              } from '../modules/combine_gvcfs' 
include { CONVERT_INTERVALS                                             } from '../modules/convert_intervals' 
include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 
include { GENOTYPE_POSTERIORS                                              } from '../modules/genotype_posteriors' 



workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed

    main: 

    // check that params.interval_size is not less than 100k bases
    if ( params.interval_size.toFloat() < 1E5 ){
        println ("\n*** ERROR: 'params.interval_size' must be >= 1E5 (100,000 bases) ***\n")
        error()
    }

    /* 
        Read QC
    */

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
        .set { ch_interval_list }

    // call variants for single samples
    CALL_VARIANTS (
        ch_sample_bam,
        ch_genome_indexed
    )

    // combine GVCFs into a single element
    CALL_VARIANTS.out.gvcf
        .map { sample, gvcf, tbi -> [ gvcf, tbi ] }
        .collect()
        .set { ch_gvcfs }

    // combine GVCFs into one file
    COMBINE_GVCFS (
        ch_gvcfs,
        ch_genome_indexed
    )

    // calculate genotype posteriors over each genomic interval
    GENOTYPE_POSTERIORS (
        COMBINE_GVCFS.out.gvcf,
        ch_interval_list
    )

    // call genotypes at variant sites

    emit: 
    CALL_VARIANTS.out


}