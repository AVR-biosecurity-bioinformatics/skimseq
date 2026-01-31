/*
    Filter .vcf files from GATK
*/

//// import modules
include { CALC_CHUNK_DP                                } from '../modules/calc_chunk_dp'
include { MERGE_CHUNK_DP                               } from '../modules/merge_chunk_dp'
include { FILTER_VCF_SITES                             } from '../modules/filter_vcf_sites'
include { EXTRACT_VCF_SITES                            } from '../modules/extract_vcf_sites'
include { MERGE_VCFS as MERGE_FILTERED_VCFS            } from '../modules/merge_vcfs'
include { MERGE_VCFS as MERGE_FILTERED_SITELISTS       } from '../modules/merge_vcfs'
include { VCF_STATS                                    } from '../modules/vcf_stats'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_VARIANTS {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_mask_bed_vcf
    ch_sample_names

    main: 

    // Calculate missing data and variant DP histogram for each chunk
    CALC_CHUNK_DP (
        ch_vcfs
    )

    // Merge output into single 
    CALC_CHUNK_DP.out.chunk_dp
            .map { interval_hash, interval_bed, bed_tbi, dphist -> dphist }
            .collect()
            .set { ch_chunk_dp }

    // Merge missing data and DP histogram from all chunks
    MERGE_CHUNK_DP (
        ch_chunk_dp
    )

    // For each input VCF, combine with type to make a copy for each variant type, then run FILTER_VCF_SITES on each
    def types = ['snp', 'indel']
    if( params.output_invariant ) {
    types << 'invariant'
    }

    FILTER_VCF_SITES (
        ch_vcfs.combine( channel.of(*types) ),
	    ch_mask_bed_vcf,
        MERGE_CHUNK_DP.out.dp_hist
    )

    // Use counts file to remove those chunks which contain no variants
    FILTER_VCF_SITES.out.vcf
        .map { type, interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> tuple(type, interval_hash, interval_bed, bed_tbi, vcf, tbi) }
        .set { ch_vcfs_nonempty }

    // Merge variant types back together, by chunk
    ch_vcfs_nonempty
        .map { type, interval_hash, interval_bed, bed_tbi, vcf, tbi ->
            tuple("${interval_hash}_filtered", vcf, tbi) // Adding string '_filtered' to ihas ensure's its present in from MERGE_FILTERED_VCF output filename
        }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_FILTERED_VCFS (
        ch_vcf_to_merge
    )
  
   // Output channel of interval_hash, interval_bed, vcf, tbi, sitesvcf, sitestbi
    ch_vcfs   
        .join ( MERGE_FILTERED_VCFS.out.vcf 
                    .map { interval_hash_filtered, vcf, tbi ->
                    def interval_hash = interval_hash_filtered.replaceFirst(/_filtered$/, '') // remove '_filtered' string from ihash for join
                    tuple(interval_hash, vcf, tbi)
                }, by: 0 )
        .set { ch_vcf_filtered }

    // Output channels of just the merged sitelists
    
    // Create a channel of all 3 variant types + all together for merging
    ch_vcfs_nonempty.map { type, interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple(type, vcf, tbi) }
        .concat(ch_vcfs_nonempty.map { type, interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple('combined', vcf, tbi) })
        .groupTuple(by: 0)
       .set { ch_sitelists_to_merge }


    // Group all filtered sitelists by variant type and merge
    MERGE_FILTERED_SITELISTS (
        ch_sitelists_to_merge
    )
   
    // Extract merged variant type vcfs into convenient channels
    MERGE_FILTERED_SITELISTS.out.vcf.filter{ it[0]=='combined' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_all_sitelist }
    MERGE_FILTERED_SITELISTS.out.vcf.filter{ it[0]=='snp' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_snp_sitelist }
    MERGE_FILTERED_SITELISTS.out.vcf.filter{ it[0]=='indel' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_indel_sitelist }
    MERGE_FILTERED_SITELISTS.out.vcf.filter{ it[0]=='invariant' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_inv_sitelist }

    // QC plots for sites and genotypes
    PLOT_VCF_FILTERS (
        FILTER_VCF_SITES.out.hist.collect(),
        FILTER_VCF_SITES.out.summary.collect()
    )

    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_combined = ch_vcf_filtered
    all_sitelist = ch_all_sitelist
    snp_sitelist = ch_snp_sitelist
    indel_sitelist = ch_indel_sitelist
    inv_sitelist = ch_inv_sitelist
    //reports = VCF_STATS.out.vcfstats

}