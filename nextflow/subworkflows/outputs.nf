/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle/create_beagle' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination/plot_ordination' 
include { PLOT_PCA                                               } from '../modules/plot_pca/plot_pca' 
include { PLOT_TREE                                              } from '../modules/plot_tree/plot_tree' 
include { PLINK_IMPORT                                           } from '../modules/plink_import/plink_import' 
include { PLINK_PCA                                              } from '../modules/plink_pca/plink_pca' 
include { PLINK_REL                                              } from '../modules/plink_rel/plink_rel' 
include { PLINK_KING                                             } from '../modules/plink_king/plink_king' 
include { PLINK_DIST                                             } from '../modules/plink_dist/plink_dist' 

workflow OUTPUTS {

    take:
    ch_final_bcfs
    ch_genome_indexed
    ch_sample_pop

    main: 

    /* 
        Create outputs
    */


    // Create beagle GL file
    ch_beagle_gl = channel.empty()
    if (params.output_beagle_gl) {
        CREATE_BEAGLE_GL (
            ch_final_bcfs,
            ch_genome_indexed,
            false
        )
        ch_beagle_gl = CREATE_BEAGLE_GL.out.beagle
    }

    // Import PLINK file
    PLINK_IMPORT (
        ch_final_bcfs
    )

    // Run PCA on plink bed
    PLINK_PCA (
        PLINK_IMPORT.out.plink
    )   

    // Create relationship matrix from plink bed
    PLINK_REL (
        PLINK_IMPORT.out.plink
    )   
    
    // Create KING relationship matrix from plink bed
    PLINK_KING (
        PLINK_IMPORT.out.plink
    )   

    // Create PLINK IBS dist matrix matrix from plink bed
    PLINK_DIST (
        PLINK_IMPORT.out.plink
    )   

    // Turn ch_sample_pop tuples into a 2‑col TSV 'popmap' file
    ch_sample_pop
        .map { sample, pop -> "$sample\t$pop" }
        .collectFile(name: 'sample_pop.tsv', newLine: true)
        .first()
        .set { ch_popmap }

    // create ordination plot from distance matrices
    PLOT_ORDINATION (
        PLINK_DIST.out.mat,
        ch_popmap,
        false
    )

    // create PCA plot from PLINK outputs
    PLOT_PCA (
        PLINK_PCA.out.pca,
        ch_popmap
    )

    // Create NJ tree from distance matrix
    PLOT_TREE (
        PLINK_DIST.out.mat,
        ch_popmap
    )



    emit:
    beagle_gl        = ch_beagle_gl
    plink            = PLINK_IMPORT.out.plink
    pca              = PLINK_PCA.out.pca
    relationship     = PLINK_REL.out.rel
    king             = PLINK_KING.out.king
    distance         = PLINK_DIST.out.mat
    ordination_plot  = PLOT_ORDINATION.out.plots
    pca_plot         = PLOT_PCA.out.plots
    tree_plot        = PLOT_TREE.out.plots
    newick_tree      = PLOT_TREE.out.newick_tree
    popmap           = ch_popmap

}