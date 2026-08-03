/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle/create_beagle' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination/plot_ordination' 
include { PLOT_PCA                                               } from '../modules/plot_pca/plot_pca' 
include { PLOT_TREE                                              } from '../modules/plot_tree/plot_tree' 
include { CONCAT_VCFS as CONCAT_FINAL                            } from '../modules/concat_vcfs/concat_vcfs'
include { SPLIT_VCF_BY_TYPE                                      } from '../modules/split_vcf_by_type/split_vcf_by_type'
include { PLINK_IMPORT                                           } from '../modules/plink_import/plink_import' 
include { PLINK_PCA                                              } from '../modules/plink_pca/plink_pca' 
include { PLINK_REL                                              } from '../modules/plink_rel/plink_rel' 
include { PLINK_KING                                             } from '../modules/plink_king/plink_king' 
include { PLINK_DIST                                             } from '../modules/plink_dist/plink_dist' 

workflow OUTPUTS {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_sample_pop

    main: 

    /* 
        Create outputs
    */

    // First split chunked vcfs by type
    SPLIT_VCF_BY_TYPE(
        ch_vcfs.map { interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple(interval_hash, vcf, tbi) }
    )

    // Build merge input channels from the named emits
    def ch_merge_inputs = SPLIT_VCF_BY_TYPE.out.snp_vcf
        .map { interval_hash, vcf, tbi -> tuple('snp', vcf, tbi) }

    if( params.output_indel ) {
        ch_merge_inputs = ch_merge_inputs.mix(
            SPLIT_VCF_BY_TYPE.out.indel_vcf
                .map { interval_hash, vcf, tbi -> tuple('indel', vcf, tbi) }
        )
    }

    if( params.output_invariant ) {
        ch_merge_inputs = ch_merge_inputs.mix(
            SPLIT_VCF_BY_TYPE.out.invariant_vcf
                .map { interval_hash, vcf, tbi -> tuple('invariant', vcf, tbi) }
        )
    }

    // Keep the combined merge from the original chunk VCFs
    ch_merge_inputs = ch_merge_inputs.mix(
        ch_vcfs.map { interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple('combined', vcf, tbi) }
    )

    // Group all chunked vcfs by variant type and merge
    ch_merge_inputs
        .groupTuple(by: 0)
        .set { ch_filtered_vcfs_to_merge }

    // Group all filtered sitelists by variant type and merge
    CONCAT_FINAL (
        ch_filtered_vcfs_to_merge
    )
   
    // Extract merged variant type vcfs into convenient channels
    CONCAT_FINAL.out.vcf.filter{ it[0]=='combined' }.first().set { ch_final_all }
    CONCAT_FINAL.out.vcf.filter{ it[0]=='snp' }.first().set { ch_final_snp }
    CONCAT_FINAL.out.vcf.filter{ it[0]=='indel' }.first().set { ch_final_indel }
    CONCAT_FINAL.out.vcf.filter{ it[0]=='invariant' }.first().set { ch_final_inv }

    /* 
        Create outputs
    */


    // Create channel containing filtered VCF along with seperate SNP and INDEL vcf
    ch_final_all
        .mix(ch_final_snp, ch_final_indel)
        .set{ ch_final_vcfs }

    // Create beagle GL file
    ch_beagle_gl = Channel.empty()
    if (params.output_beagle_gl) {
        CREATE_BEAGLE_GL (
            ch_final_vcfs,
            ch_genome_indexed,
            false
        )
        ch_beagle_gl = CREATE_BEAGLE_GL.out.beagle
    }

    // Import PLINK file
    PLINK_IMPORT (
        ch_final_vcfs
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
        .map { s,p -> "$s\t$p" }
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

    // Create NJ tree
    PLOT_TREE (
        PLINK_DIST.out.mat,
        ch_popmap
    )

    emit:
    final_vcf_all    = ch_final_all.map{ name, vcf, tbi -> tuple( vcf, tbi)}
    final_vcf        = CONCAT_FINAL.out.vcf
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