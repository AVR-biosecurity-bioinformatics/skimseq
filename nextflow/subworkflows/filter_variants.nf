/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_INDELS                                                 } from '../modules/filter_indels'
include { FILTER_SNPS                                                 } from '../modules/filter_snps'
include { MERGE_FILTERED                                                 } from '../modules/merge_filtered'



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

    MERGE_FILTERED (
        FILTER_SNPS.out.vcf,
        FILTER_INDELS.out.vcf
    )



    emit: 
    FILTER_SNPS.out.vcf


}