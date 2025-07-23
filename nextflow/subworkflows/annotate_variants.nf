/*
    annotate .vcf files for filtering
*/

//// import modules
include { ANNOTATE_VCF                                            } from '../modules/annotate_vcf'

workflow FILTER_VARIANTS {

    take:
    ch_vcf

    main: 

    // annotate VCF
    ANNOTATE_VCF (
        ch_vcf
    )
           
    emit: 
    annotated_vcf = ANNOTATE_VCF.out.vcf

}