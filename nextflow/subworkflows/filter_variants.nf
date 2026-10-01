/*
    Filter .vcf files 
*/

//// import modules
include { CALC_CHUNK_DP                                } from '../modules/calc_chunk_dp/calc_chunk_dp'
include { MERGE_CHUNK_DP                               } from '../modules/merge_chunk_dp/merge_chunk_dp'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing/merge_chunk_missing'
include { FILTER_VCF                                   } from '../modules/filter_vcf/filter_vcf'
include { PLOT_VARIANT_FILTERS                         } from '../modules/plot_variant_filters/plot_variant_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters/plot_sample_filters'

workflow FILTER_VARIANTS {

    take:
    ch_bcfs
    ch_mask_bed
    ch_popmap

    main: 
   
    /*
        Calculate depth and per-sample missing data filters
    */

    // Calculate missing data and variant DP histogram for each chunk
    CALC_CHUNK_DP(
        ch_bcfs
    )

    // Merge all chunk DP histograms together
    MERGE_CHUNK_DP(
        CALC_CHUNK_DP.out.chunk_dp.map { _interval_hash, _interval_bed, _bed_tbi, dphist -> dphist }.collect(),
        params.vcf_dp_percentile_lower,
        params.vcf_dp_percentile_upper
    )

    MERGE_CHUNK_DP.out.dp_bounds
        .map { f ->
            def lines = f.readLines()
            def hdr = lines[0].split('\t')
            def row = lines[1].split('\t')
            def m = [hdr, row].transpose().collectEntries { k, v -> [(k): v] }
            tuple(m.DPlower as Integer, m.DPupper as Integer)
        }
        .set { ch_dp_bounds }

    // Merge per-sample missing data from all chunks into a single table
    MERGE_CHUNK_MISSING(
        CALC_CHUNK_DP.out.chunk_missing.map { _interval_hash, _interval_bed, _bed_tbi, missing -> missing }.collect()
    )

    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS(
        MERGE_CHUNK_MISSING.out.missing_summary
    )

    /*
        Filter VCF
    */

    // Global site filters
    FILTER_VCF(
        ch_bcfs.combine(ch_dp_bounds),
        ch_mask_bed,
        ch_popmap.first(),
        MERGE_CHUNK_MISSING.out.missing_summary
    )

    // Remove chunks which contain no variants after filtering
    FILTER_VCF.out.bcf
        .map { interval_hash, interval_bed, bed_tbi, bcf, csi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(interval_hash, interval_bed, bed_tbi, bcf, csi, n)
        }
        .filter { _interval_hash, _interval_bed, _bed_tbi, _bcf, _csi, n -> n > 0 }
        .map { interval_hash, interval_bed, bed_tbi, bcf, csi, _n ->
            tuple(interval_hash, interval_bed, bed_tbi, bcf, csi)
        }
        .set { ch_filtered_bcf }

    // Create list of samples surviving filtering
    FILTER_VCF.out.samples_to_keep.first()
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }


    // QC plots for site histograms
    PLOT_VARIANT_FILTERS (
        FILTER_VCF.out.metrics.map { _interval_hash, _interval_bed, _bed_tbi, tsv -> tsv }.collect(),
        "site_filters"
    )

    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_bcf = ch_filtered_bcf
    sample_names_filt = ch_sample_names_filt
    sample_filter_plots = PLOT_SAMPLE_FILTERS.out.plots
    sample_missing_tsv = PLOT_SAMPLE_FILTERS.out.sample_missing_tsv
    site_filter_plots = PLOT_VARIANT_FILTERS.out.plots
}