/*
    annotate .vcf files for filtering
*/

//// import modules
include { ANNOTATE_VCF                                            } from '../modules/annotate_vcf'

workflow ANNOTATE_VARIANTS {

    take:
    ch_vcf

    main: 

    // TODO: Process to estimate AF from GL's will go here: issue #103


    // annotate VCF
    ANNOTATE_VCF (
        ch_vcf
    )
           
    emit: 
    vcf = ANNOTATE_VCF.out.vcf

}