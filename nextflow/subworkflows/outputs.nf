/*
    Create outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
include { VCF2DIST                                               } from '../modules/vcf2dist' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination' 
include { PLOT_TREE                                              } from '../modules/plot_tree' 
include { MERGE_VCFS as MERGE_FINAL                              } from '../modules/merge_vcfs'
include { SPLIT_VCF_BY_TYPE                                      } from '../modules/split_vcf_by_type'

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
    MERGE_FINAL (
        ch_filtered_vcfs_to_merge
    )
   
    // Extract merged variant type vcfs into convenient channels
    MERGE_FINAL.out.vcf.filter{ it[0]=='combined' }.first().set { ch_final_all }
    MERGE_FINAL.out.vcf.filter{ it[0]=='snp' }.first().set { ch_final_snp }
    MERGE_FINAL.out.vcf.filter{ it[0]=='indel' }.first().set { ch_final_indel }
    MERGE_FINAL.out.vcf.filter{ it[0]=='invariant' }.first().set { ch_final_inv }

    /* 
        Create outputs
    */


    // Create channel containing filtered VCF along with seperate SNP and INDEL vcf
    ch_final_all
        .mix(ch_final_snp, ch_final_indel)
        .set{ ch_final_vcfs }

    // Create beagle GL file
    ch_beagle_gl_out = Channel.empty()
    if (params.output_beagle_gl) {
        CREATE_BEAGLE_GL (
            ch_final_vcfs,
            ch_genome_indexed,
            false
        )
        ch_beagle_gl_out = CREATE_BEAGLE_GL.out.beagle
    }

    // Create distance matrices from VCFs
    VCF2DIST (
        ch_final_vcfs
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
    vcf = ch_final_all.map{ name, vcf, tbi -> tuple( vcf, tbi)}
}