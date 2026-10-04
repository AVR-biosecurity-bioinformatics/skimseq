/*
    Filter .vcf files 
*/

//// import modules
include { MERGE_RIKER_COVERAGE                         } from '../modules/merge_riker_coverage/merge_riker_coverage'
include { FILTER_VCF                                   } from '../modules/filter_vcf/filter_vcf'
include { PLOT_VARIANT_FILTERS                         } from '../modules/plot_variant_filters/plot_variant_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters/plot_sample_filters'
include { CONCAT_VCFS as CONCAT_FINAL                  } from '../modules/concat_vcfs/concat_vcfs'
include { VCF_STATS as VCF_STATS_FILTERED              } from '../modules/vcf_stats/vcf_stats'

workflow FILTER_VARIANTS {

    take:
    ch_bcfs
    ch_mask_bed
    ch_popmap
    ch_wgs_coverage
    ch_genome_indexed

    main: 
   
    /*
        Calculate depth and per-sample missing data filters
    */

    // Merge per-sample missing data from all chunks into a single table
    ch_coverage_files = ch_wgs_coverage
        .map { _sample, stats -> stats }
        .collect()

    MERGE_RIKER_COVERAGE(
        ch_coverage_files
        )

    ch_missing_summary = MERGE_RIKER_COVERAGE.out.missing_summary
    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS(
        ch_missing_summary
    )

    /*
        Filter VCF
    */

    // Global site filters
    FILTER_VCF(
        ch_bcfs,
        ch_mask_bed,
        ch_popmap.first(),
        ch_missing_summary
    )

    // Create list of samples surviving filtering
    FILTER_VCF.out.samples_to_keep.first()
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }

    ch_filter_hists = FILTER_VCF.out.filter_hist
        .map { _interval_hash, histogram -> histogram }
        .collect()

    // QC plots for site histograms
    PLOT_VARIANT_FILTERS(ch_filter_hists)

    // Build merge input channels from the named emits
    def ch_merge_inputs = FILTER_VCF.out.snp_bcf
        .map { _interval_hash, _interval_bed, _bed_tbi, bcf, csi -> tuple('snp', bcf, csi) }

    if( params.output_indel ) {
        ch_merge_inputs = ch_merge_inputs.mix(
            FILTER_VCF.out.indel_bcf
                .map { _interval_hash, _interval_bed, _bed_tbi, bcf, csi -> tuple('indel', bcf, csi) }
        )
    }

    if( params.output_invariant ) {
        ch_merge_inputs = ch_merge_inputs.mix(
            FILTER_VCF.out.invariant_bcf
                .map { _interval_hash, _interval_bed, _bed_tbi, bcf, csi -> tuple('invariant', bcf, csi) }
        )
    }

    // Keep the combined merge from the original chunk VCFs
    ch_merge_inputs = ch_merge_inputs.mix(
            FILTER_VCF.out.all_bcf
                .map { _interval_hash, _interval_bed, _bed_tbi, bcf, csi -> tuple('combined', bcf, csi) }
        )


    // Group all chunked vcfs by variant type and merge
    ch_merge_inputs
        .groupTuple(by: 0)
        .set { ch_filtered_vcfs_to_merge }

    // Group all filtered sitelists by variant type and merge
    CONCAT_FINAL (
        ch_filtered_vcfs_to_merge
    )

    // Calculate VCF statistics on the final file
    VCF_STATS_FILTERED (
        CONCAT_FINAL.out.vcf.filter { record -> record[0] == 'combined' }.map{ _name, bcf, csi -> tuple( bcf, csi)},
        ch_genome_indexed
    )

    // Subset the merged vcf channels to each variant type for emission
    emit:
    final_bcf = CONCAT_FINAL.out.vcf
    sample_names_filt = ch_sample_names_filt
    sample_filter_plots = PLOT_SAMPLE_FILTERS.out.plots
    missing_summary = ch_missing_summary
    site_filter_plots = PLOT_VARIANT_FILTERS.out.plots
    vcf_stats = VCF_STATS_FILTERED.out.vcfstats
}