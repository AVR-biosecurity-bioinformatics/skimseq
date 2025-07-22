/*
    Low coverage outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
include { CREATE_BEAGLE as CREATE_BEAGLE_GP                      } from '../modules/create_beagle' 
include { CREATE_PSEUDOHAP                                       } from '../modules/create_pseudohap' 
include { VCF2DIST                                               } from '../modules/vcf2dist' 

workflow LOWCOV_OUTPUTS {

    take:
    ch_filtered_vcf
    ch_genome_indexed

    main: 

    /* 
        Probablistic genotyping
    */

    // Create beagle GL file
    CREATE_BEAGLE_GL (
        ch_filtered_vcf,
        ch_genome_indexed,
        false
    )

    // Create beagle GP file
    CREATE_BEAGLE_GP (
        ch_filtered_vcf,
        ch_genome_indexed,
        true
    )

    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
        ch_filtered_vcf,
        ch_genome_indexed
    )

    // Create distance matrices
    ch_filtered_vcf
        .mix(CREATE_PSEUDOHAP.out.vcf)
        .set{ ch_vcfs }

    ch_vcfs.view()

    VCF2DIST (
        ch_vcfs
    )

    // TODO - create ordination plot
    // This should take in a distance matrix or covariant matrix
    // Use any population labels specified in the input file for labelling and colours


    // TODO - Create NJ tree
    // This should take in just a distance matrix
    // Use any population labels specified in the input file for labelling and colours

    emit: 
    beagle_gl = CREATE_BEAGLE_GL.out.beagle
    beagle_gp = CREATE_BEAGLE_GP.out.beagle
    pseudo = CREATE_PSEUDOHAP.out.tsv


}