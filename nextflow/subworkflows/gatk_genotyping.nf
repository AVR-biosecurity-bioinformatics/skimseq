/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                             } from '../modules/call_variants'


workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed

    main: 

    /* 
        Read QC
    */

    CALL_VARIANTS (
        ch_sample_bam,
        ch_genome_indexed
    )

    

    emit: 
    CALL_VARIANTS.out


}