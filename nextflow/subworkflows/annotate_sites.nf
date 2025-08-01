/*
    annotate .vcf files for filtering
*/

//// import modules
include { ANNOTATE_VCF                                            } from '../modules/annotate_vcf'

workflow ANNOTATE_SITES {

    take:
    ch_vcf

    main: 

    // TODO: Process to estimate AF from GL's will go here: issue #103

    // TODO: Maybe worth making the annotation files in separate processes, with annotate vcfs just adding them

    // annotate VCF
    ANNOTATE_VCF (
        ch_vcf
    )
           
    emit: 
    vcf = ANNOTATE_VCF.out.vcf

}