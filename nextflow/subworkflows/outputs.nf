/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
include { CREATE_BEAGLE as CREATE_BEAGLE_GP                      } from '../modules/create_beagle' 
include { CREATE_PSEUDOHAP                                       } from '../modules/create_pseudohap' 
include { VCF2DIST                                               } from '../modules/vcf2dist' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination' 
include { PLOT_TREE                                              } from '../modules/plot_tree' 

workflow OUTPUTS {

    take:
    ch_filtered_vcf
    ch_genome_indexed
    ch_sample_pop

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

    VCF2DIST (
        ch_vcfs
    )

    // create ordination plot
    // This should take in a distance matrix or covariant matrix
    // TODO: Use any population labels specified in the input file for labelling and colours
    PLOT_ORDINATION (
        VCF2DIST.out.mat,
        ch_sample_pop,
        false
    )

    // Create NJ tree
    // TODO: Use any population labels specified in the input file for labelling and colours
    //PLOT_TREE (
    //    VCF2DIST.out.mat
    //)

    emit: 
    beagle_gl = CREATE_BEAGLE_GL.out.beagle
    beagle_gp = CREATE_BEAGLE_GP.out.beagle
    pseudo = CREATE_PSEUDOHAP.out.tsv


}