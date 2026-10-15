nextflow.enable.types = true

/*
 * Create downstream analysis outputs.
 */

include { CREATE_BEAGLE as CREATE_BEAGLE_GL } from '../modules/create_beagle/create_beagle'
include { PLOT_ORDINATION                   } from '../modules/plot_ordination/plot_ordination'
include { PLOT_PCA                          } from '../modules/plot_pca/plot_pca'
include { PLOT_TREE                         } from '../modules/plot_tree/plot_tree'
include { PLINK_IMPORT                      } from '../modules/plink_import/plink_import'
include { PLINK_PCA                         } from '../modules/plink_pca/plink_pca'
include { PLINK_REL                         } from '../modules/plink_rel/plink_rel'
include { PLINK_KING                        } from '../modules/plink_king/plink_king'
include { PLINK_DIST                        } from '../modules/plink_dist/plink_dist'

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
     * Create Beagle genotype-likelihood output.
     */

    beagle_gl = channel.empty()

    if (params.output_beagle_gl) {
        beagle_gl = CREATE_BEAGLE_GL(
            snp_bcf,
            false
        )
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

    pca = PLINK_PCA(
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

    pca_plot = PLOT_PCA(
        pca,
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
     * Outputs.
     */

    emit:
    beagle_gl = beagle_gl
    plink: Channel<PlinkDataset> = plink
    pca = pca
    relationship = relationship
    king = king
    distance = distance
    ordination_plot = ordination_plot
    pca_plot = pca_plot
    tree = tree
}