nextflow.enable.types = true

/*
 * Create downstream analysis outputs.
 */

include { PLOT_ORDINATION                   } from '../modules/plot_ordination/plot_ordination'
include { PLOT_PCA as PLOT_PLINK_PCA        } from '../modules/plot_pca/plot_pca'
include { PLOT_PCA as PLOT_PCAONE_PCA       } from '../modules/plot_pca/plot_pca'
include { PLOT_TREE                         } from '../modules/plot_tree/plot_tree'
include { PLINK_IMPORT                      } from '../modules/plink_import/plink_import'
include { PLINK_PCA                         } from '../modules/plink_pca/plink_pca'
include { PLINK_REL                         } from '../modules/plink_rel/plink_rel'
include { PLINK_KING                        } from '../modules/plink_king/plink_king'
include { PLINK_DIST                        } from '../modules/plink_dist/plink_dist'
include { CREATE_BEAGLE                     } from '../modules/create_beagle/create_beagle'
include { PCAONE_RELATEDNESS                } from '../modules/pcaone_relatedness/pcaone_relatedness'
include { PCAONE_PCA                        } from '../modules/pcaone_pca/pcaone_pca'

// Import record types
include {
    VariantFile
    PlinkDataset
} from '../types'

workflow OUTPUTS {

    take:
    final_bcfs: Channel<VariantFile>
    popmap: Value<Path>

    main:

    /*
     * Select the final SNP dataset.
     */

    snp_bcf = final_bcfs
        .filter { result ->
            result.outname == 'snp'
        }


    /*
     * Import filtered SNPs into PLINK format.
     */

    plink = PLINK_IMPORT(
        snp_bcf
    )


    /*
     * PLINK analyses.
     */

    plink_pca = PLINK_PCA(
        plink
    )

    relationship = PLINK_REL(
        plink
    )

    king = PLINK_KING(
        plink
    )

    distance = PLINK_DIST(
        plink
    )


    /*
     * Plot distance-matrix ordination.
     */

    ordination_plot = PLOT_ORDINATION(
        distance,
        popmap,
        false
    )


    /*
     * Plot PCA.
     */

    pca_plot_plink = PLOT_PLINK_PCA(
        plink_pca,
        popmap
    )


    /*
     * Generate neighbour-joining tree.
     */

    tree = PLOT_TREE(
        distance,
        popmap
    )


    /*
     * PCAone analyses
     */

    beagle = CREATE_BEAGLE(
        snp_bcf,
        false
    )

    pcaone_pca = PCAONE_PCA(
        beagle
    )

    pca_plot_pcaone = PLOT_PCAONE_PCA(
        pcaone_pca,
        popmap
    )

    pcaone_rel = PCAONE_RELATEDNESS(
        plink,
        params.relatedness_n_pcs
    )


    /*
     * Outputs.
     */

    emit:
    beagle_gl = beagle
    pcaone_rel = pcaone_rel
    pcaone_pca = pcaone_pca
    plink: Channel<PlinkDataset> = plink
    plink_pca = plink_pca
    relationship = relationship
    king = king
    distance = distance
    ordination_plot = ordination_plot
    pca_plot_plink = pca_plot_plink
    pca_plot_pcaone = pca_plot_pcaone
    tree = tree
}