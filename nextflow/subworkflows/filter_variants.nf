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
    ch_sample_pop

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
    def variant_types = ['snp', 'indel']
    if( params.output_invariant ) {
    variant_types << 'invariant'
    }

   channel.of(*variant_types) 
	.combine(ch_vcfs)
	.set { ch_vcf_types }

    // Global sites filters - TODO: this adds tags to the vcf but does not filter
    FILTER_VCF_SITES (
        ch_vcf_types,
	    ch_mask_bed_vcf,
        MERGE_CHUNK_DP.out.dp_hist
    )

    // Create seperate channels that have all sample names grouped by pop

    // Filter on a per-populations basis, outputting a sites only vcf of just passing sites
    //FILTER_VCF_POP (
    //    ch_vcf_types,
    //    ch_sample_pop
    //)

    // intersect per-pop sitelists, keeping only those failing in >n samples
    // Then intersect with global filter


    // Create site histograms - uses the tages from the soft filtered vcf
    //CREATE_FILTER_HISTS(
    //    
    //)

    // QC plots for site histograms
    //PLOT_VCF_FILTERS (
    //    FILTER_VCF_SITES.out.hist.collect(),
    //    FILTER_VCF_SITES.out.summary.collect(),
    //    "site_filters"
    //)


    // Use counts file to remove those chunks which contain no variants
    FILTER_VCF_SITES.out.vcf
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> tuple(variant_type, interval_hash, vcf, tbi) }
        .set { ch_vcfs_nonempty }


   // Output channel of variant_type, interval_hash, interval_bed, vcf, tbi, sitesvcf, sitestbi
    ch_vcf_types   
        .join(ch_vcfs_nonempty, by: [0,1] )
        .set { ch_filtered_sitelist }

    // Output channels of just the merged sitelists
    
    // Create a channel of all 3 variant types + all together for merging
    ch_vcfs_nonempty.map { variant_type, interval_hash, vcf, tbi -> tuple(variant_type, vcf, tbi) }
        .concat(ch_vcfs_nonempty.map { variant_type, interval_hash, vcf, tbi -> tuple('combined', vcf, tbi) })
        .groupTuple(by: 0)
        .set { ch_sitelists_to_merge }


    // Group all filtered sitelists by variant type and merge
    MERGE_FILTERED_SITELISTS (
        ch_sitelists_to_merge
    )
   
    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_sitelist = ch_filtered_sitelist

}