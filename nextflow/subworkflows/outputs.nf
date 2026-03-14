/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
include { VCF2DIST                                               } from '../modules/vcf2dist' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination' 
include { PLOT_TREE                                              } from '../modules/plot_tree' 
include { MERGE_VCFS as MERGE_FINAL                              } from '../modules/merge_vcfs'

workflow OUTPUTS {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_sample_pop

    main: 

    /* 
        Create outputs
    */

    // TODO: need to split vcf into variant types here

    // Define the three variant types with optional inclusion of indel and invariant
    def variant_types = ['snp']
    if( params.output_indel )      variant_types << 'indel'
    if( params.output_invariant )  variant_types << 'invariant'


    ch_vcfs
        .combine(Channel.of(*variant_types))
        .set { ch_vcf_types }

    // Create a channel of all 3 variant types + all together for merging
     ch_vcfs
        .combine(Channel.of(*variant_types))
        ch_vcfs.map { interval_hash, interval_bed, bed_tbi, vcf, tbi, type -> tuple(type, vcf, tbi) }
        .concat(ch_vcfs.map { interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple('combined', vcf, tbi) })
        .groupTuple(by: 0)
        .set { ch_filtered_vcfs_to_merge }

    // Group all filtered sitelists by variant type and merge
    MERGE_FINAL (
        ch_filtered_vcfs_to_merge
    )
   
    // Extract merged variant type vcfs into convenient channels
    MERGE_FINAL.out.vcf.filter{ it[0]=='combined' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_final_all }
    MERGE_FINAL.out.vcf.filter{ it[0]=='snp' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_final_snp }
    MERGE_FINAL.out.vcf.filter{ it[0]=='indel' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_final_indel }
    MERGE_FINAL.out.vcf.filter{ it[0]=='invariant' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_final_inv }

    /* 
        Create outputs
    */


    // Create channel containing filtered VCF along with seperate SNP and INDEL vcf
    ch_final_all
        //.mix(ch_final_snp, ch_final_indel)
        .set{ ch_final_vcfs }

    // Create beagle GL file
    def ch_beagle_gl_out = Channel.empty()
    if (params.output_beagle_gl) {
        CREATE_BEAGLE_GL (
            ch_final_vcfs,
            ch_genome_indexed,
            false
        )
        ch_beagle_gl_out = CREATE_BEAGLE_GL.out.beagle
    }

    // Create updated channel for distance matrices
    ch_final_vcfs
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
    emit:
    vcf = ch_final_all
}