/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_SNPS                                                 } from '../modules/filter_snps'



workflow FILTER_VARIANTS {

    take:
    ch_vcf

    main: 

    FILTER_SNPS (
        ch_vcf
    )
    

    emit: 
    FILTER_SNPS.out.vcf


}