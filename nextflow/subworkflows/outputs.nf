/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
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
        //.mix(ch_filtered_snps, ch_filtered_indels)
        .set{ ch_vcfs }

    // Create beagle GL file
    def ch_beagle_gl_out = Channel.empty()
    if (params.output_beagle_gl) {
        CREATE_BEAGLE_GL (
            ch_vcfs,
            ch_genome_indexed,
            false
        )
        ch_beagle_gl_out = CREATE_BEAGLE_GL.out.beagle
    }

    // Create updated channel for distance matrices
    ch_vcfs
        .set{ ch_vcfs_for_dist }

    // Create distance matrices from VCFs
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

    // Create NJ tree
    PLOT_TREE (
        VCF2DIST.out.mat,
        ch_popmap
    )

}