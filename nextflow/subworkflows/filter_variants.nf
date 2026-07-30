/*
    Filter .vcf files from GATK
*/

//// import modules
include { EXTRACT_VCF_SITES                            } from '../modules/extract_vcf_sites/extract_vcf_sites'
include { COUNT_VCF_RECORDS                            } from '../modules/count_vcf_records/count_vcf_records'
include { SUBSET_VCF_TO_SITES                          } from '../modules/subset_vcf_to_sites/subset_vcf_to_sites'
include { CALC_CHUNK_DP                                } from '../modules/calc_chunk_dp/calc_chunk_dp'
include { MERGE_CHUNK_DP                               } from '../modules/merge_chunk_dp/merge_chunk_dp'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing/merge_chunk_missing'
include { FILTER_VCF                                   } from '../modules/filter_vcf/filter_vcf'
include { CREATE_FILTER_HIST                           } from '../modules/create_filter_hist/create_filter_hist'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters/plot_sample_filters'

workflow FILTER_VARIANTS {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_vcf
    ch_sample_names
    ch_popmap

    main: 
   
    /*
        Calculate depth and per-sample missing data filters
    */

    // Calculate missing data and variant DP histogram for each chunk
    CALC_CHUNK_DP(
        ch_vcfs
    )

    // Merge all chunk DP histograms together
    MERGE_CHUNK_DP(
        CALC_CHUNK_DP.out.chunk_dp.map { interval_hash, interval_bed, bed_tbi, dphist -> dphist }.collect(),
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
        CALC_CHUNK_DP.out.chunk_missing.map { interval_hash, interval_bed, bed_tbi, missing -> missing }.collect()
    )

    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS(
        MERGE_CHUNK_MISSING.out.missing_summary,
        params.vcf_sample_max_missing
    )

    /*
        Filter VCF
    */

    // Global site filters
    FILTER_VCF(
        ch_vcfs.combine(ch_dp_bounds),
        ch_mask_bed_vcf,
        ch_popmap.first(),
        MERGE_CHUNK_MISSING.out.missing_summary
    )

    // Remove chunks which contain no variants after filtering
    FILTER_VCF.out.vcf
        .map { interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { interval_hash, interval_bed, bed_tbi, vcf, tbi, n ->
            tuple(interval_hash, interval_bed, bed_tbi, vcf, tbi)
        }
        .set { ch_filtered_vcf }

     /*
        Create site filtering QC plots
    */

    // QC plots for site histograms
    PLOT_VCF_FILTERS (
        FILTER_VCF.out.metrics.map { interval_hash, interval_bed, bed_tbi, tsv -> tsv }.collect(),
        "site_filters"
    )

     /*
        Create sitelist files for re-genotyping
    */
    
    FILTER_VCF.out.sitelist
        .map { interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple( interval_hash, interval_bed, bed_tbi, vcf, tbi, n )
        }
        .filter { interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> tuple( interval_hash, interval_bed, bed_tbi, vcf, tbi) }
        .set { ch_filtered_sites }

    FILTER_VCF.out.samples_to_keep.first()
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }
   
    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_vcf = ch_filtered_vcf
    filtered_sitelist = ch_filtered_sites
    sample_names_filt = ch_sample_names_filt
}