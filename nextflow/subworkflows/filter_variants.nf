/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_INDELS                                                 } from '../modules/filter_indels'
include { FILTER_SNPS                                                 } from '../modules/filter_snps'



workflow FILTER_VARIANTS {

    take:
    ch_vcf

    main: 

    FILTER_SNPS (
        ch_vcf
    )
    
    FILTER_INDELS (
        ch_vcf
    )



    emit: 
    FILTER_SNPS.out.vcf


}