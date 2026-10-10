nextflow.enable.types = true

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

// Import record types
include {
    VariantChunk
    VariantFile
    PreparedGenome
    CramStats
} from '../types'


workflow FILTER_VARIANTS {

    take:
    bcfs: Channel<VariantChunk>
    mask_bed: Value<Path>
    popmap: Value<Path>
    wgs_coverage: Channel<CramStats>
    genome: Value<PreparedGenome>

    main: 
   
    /*
        Calculate depth and per-sample missing data filters
    */

    // Merge per-sample missing data from all chunks into a single table
    coverage_files = wgs_coverage
        .map { coverage ->
            coverage.wgs_coverage
        }
        .collect()

    missing_summary = MERGE_RIKER_COVERAGE(
        coverage_files
    )

    // QC plots for sample missing data
    sample_filter_results = PLOT_SAMPLE_FILTERS(
        missing_summary
    )
    sample_filter_plots = sample_filter_results
        .map { result -> result.plots }

    /*
        Filter VCF
    */

    // Global site filters
    filtered = FILTER_VCF(
        bcfs,
        mask_bed,
        popmap,
        missing_summary
    )

    filter_hists = filtered
        .map { chunk ->
            chunk.filter_hist
        }
        .collect()

    // QC plots for site histograms
    site_filter_plots = PLOT_VARIANT_FILTERS(
        filter_hists
    )

    // Create channel of samples surviving filtering.
    sample_names_filt = filtered
        .map { result ->
            result.samples_to_keep
        }
        .collect()
        .map { files ->
            files.find { _file ->
                true
            }
        }


    /*
    * Expand each filtered interval into the variant outputs that should
    * be concatenated.
    */

 merge_inputs = filtered
    .flatMap { result ->

        def required = [
            record(
                outname: 'snp',
                variant: result.snp_bcf,
                index: result.snp_bcf_index
            ),
            record(
                outname: 'combined',
                variant: result.all_bcf,
                index: result.all_bcf_index
            )
        ]

        def indel = params.output_indel
            ? [
                record(
                    outname: 'indel',
                    variant: result.indel_bcf,
                    index: result.indel_bcf_index
                )
            ]
            : []

        def invariant = params.output_invariant
            ? [
                record(
                    outname: 'invariant',
                    variant: result.invariant_bcf,
                    index: result.invariant_bcf_index
                )
            ]
            : []

        required + indel + invariant
    }

    /*
    * Group chunk-level variant files by output type.
    *
    * groupBy() requires keyed tuples, so temporarily convert:
    *
    *     record(outname, variant, index)
    *
    * into:
    *
    *     tuple(outname, record(variant, index))
    */

    variants_by_type = merge_inputs
        .map { result ->
            tuple(
                result.outname,
                record(
                    variant: result.variant,
                    index: result.index
                )
            )
        }
        .groupBy()

    /*
    * Convert each grouped value into the generic CONCAT_VCFS input:
    *
    *     record(
    *         outname: String,
    *         variants: Iterable<Path>,
    *         indices: Iterable<Path>
    *     )
    */

    variants_to_concat = variants_by_type
        .map { outname, files ->
            record(
                outname: outname,
                variants: files.collect { result ->
                    result.variant
                },
                indices: files.collect { result ->
                    result.index
                }
            )
        }

    /*
    * Concatenate all chunks independently for each variant type.
    */

    final_variants = CONCAT_FINAL(variants_to_concat)

    // Calculate VCF statistics on the final file
    combined_variants = final_variants
        .filter { result ->
            result.outname == 'combined'
        }
    vcf_stats = VCF_STATS_FILTERED(
        combined_variants,
        genome
    )

    // Subset the merged vcf channels to each variant type for emission
    emit:
    final_bcf: Channel<VariantFile> = final_variants
    sample_names_filt: Value<Path> = sample_names_filt
    sample_filter_plots: Value<Set<Path>> = sample_filter_plots
    missing_summary: Value<Path> = missing_summary
    site_filter_plots: Value<Set<Path>> = site_filter_plots
    vcf_stats: Channel<Path> = vcf_stats
}