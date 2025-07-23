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
    ch_filtered_merged
    ch_filtered_snps
    ch_filtered_indels
    ch_genome_indexed
    ch_sample_pop

    main: 

    /* 
        Create outputs
    */

    // Create channel containing filtered VCF along with seperate SNP and INDEL vcf
    ch_filtered_merged
        .mix(ch_filtered_snps, ch_filtered_indels)
        .set{ ch_vcfs }

    // Create beagle GL file
    CREATE_BEAGLE_GL (
        ch_vcfs,
        ch_genome_indexed,
        false
    )

    // Create beagle GP file
    CREATE_BEAGLE_GP (
        ch_vcfs,
        ch_genome_indexed,
        true
    )

    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
        ch_vcfs,
        ch_genome_indexed
    )

    // Create updated channel for distance matrices
    ch_vcfs
        .mix(CREATE_PSEUDOHAP.out.vcf)
        .set{ ch_vcfs_for_dist }

    VCF2DIST (
        ch_vcfs_for_dist
    )

    // Turn ch_sample_pop tuples into a 2‑col TSV 'popmap' file
    ch_sample_pop
        .map { s,p -> "$s\t$p" }
        .collectFile(name: 'sample_pop.tsv', newLine: true)
        .first()
        .set { ch_popmap }

    // create ordination plot from distance matrices
    PLOT_ORDINATION (
        VCF2DIST.out.mat,
        ch_popmap,
        false
    )

    // create ordination plot from covariance matrices generated via PCA
    //PLOT_ORDINATION (
    //    VCF2DIST.out.mat,
    //    ch_popmap,
    //    true
    //)

    // Create NJ tree
    PLOT_TREE (
        VCF2DIST.out.mat,
        ch_popmap
    )

    emit: 
    beagle_gl = CREATE_BEAGLE_GL.out.beagle
    beagle_gp = CREATE_BEAGLE_GP.out.beagle
    pseudo = CREATE_PSEUDOHAP.out.tsv


}