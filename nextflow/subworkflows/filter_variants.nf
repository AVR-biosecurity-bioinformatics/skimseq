/*
    Filter .vcf files from GATK
*/

//// import modules
include { EXTRACT_VCF_SITES                            } from '../modules/extract_vcf_sites'
include { COUNT_VCF_RECORDS                            } from '../modules/count_vcf_records'
include { SUBSET_VCF_TO_SITES                          } from '../modules/subset_vcf_to_sites'
include { CALC_CHUNK_DP                                } from '../modules/calc_chunk_dp'
include { MERGE_CHUNK_DP                               } from '../modules/merge_chunk_dp'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing'
include { FILTER_VCF                                   } from '../modules/filter_vcf'
include { MERGE_VCFS as MERGE_FILTERED_SITELISTS       } from '../modules/merge_vcfs'
include { CREATE_FILTER_HIST                           } from '../modules/create_filter_hist'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_VARIANTS {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_vcf
    ch_sample_names
    ch_sample_pop

    main: 

    // Define the three variant types with optional inclusion of indel and invariant
    def variant_types = ['snp']
    if( params.output_indel )      variant_types << 'indel'
    if( params.output_invariant )  variant_types << 'invariant'

    // Duplicate vcf files by variant type
    Channel.of(*variant_types)
        .combine(ch_vcfs)
        .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi ->
        tuple(vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi)
        }
        .set { ch_vcf_types }
   
    /*
        Calculate depth filters
    */


    // Calculate missing data and variant DP histogram for each chunk
    CALC_CHUNK_DP (
        ch_vcf_types
    )

    // Function to switch out dp lower and dp upper by variant type
    def dpPercForType = { String vt ->
        switch(vt) {
            case 'snp':      return [params.snp_global_dp_lower_perc,  params.snp_global_dp_upper_perc]
            case 'indel':    return [params.indel_global_dp_lower_perc,params.indel_global_dp_upper_perc]
            case 'invariant':return [params.inv_global_dp_lower_perc,  params.inv_global_dp_upper_perc]
            default:         return [null, null]
        }
    }
    // Merge chunk_dp by varaint type
    CALC_CHUNK_DP.out.chunk_dp
        .map { variant_type, interval_hash, interval_bed, bed_tbi, dphist -> tuple(variant_type, dphist ) }
        .groupTuple()
        .map { vt, files ->
        def (lo, hi) = dpPercForType(vt)
            tuple(vt, files, lo, hi)
        }
        .set { ch_chunk_dp_grouped }

    MERGE_CHUNK_DP (
        ch_chunk_dp_grouped
    )

    MERGE_CHUNK_DP.out.dp_bounds.map { vt, f ->
        def lines = f.readLines()
        def hdr = lines[0].split('\t')
        def row = lines[1].split('\t')
        def m = [hdr, row].transpose().collectEntries { k,v -> [(k): v] }
        tuple(vt, m.DPlower as Integer, m.DPupper as Integer)
    }
    .set { ch_dp_bounds }

    // Merge per-sample missing data and DP histogram from all chunks into a single table
    MERGE_CHUNK_MISSING (
         CALC_CHUNK_DP.out.chunk_missing.map { variant_type, interval_hash, interval_bed, bed_tbi, missing -> missing }.collect()
    )

    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS (
        MERGE_CHUNK_MISSING.out.missing_summary,
        params.sample_max_missing
    )

    /*
        Filter VCF
    */

    // Set up default emtpy filters
    def FILTER_DEFAULTS = [
        QUAL_THR:'NA', DPlower:'NA', PCT_LOW:'NA', DPupper:'NA', DIST_INDEL:'NA',
        EH:'NA', HWE:'NA', MAF:'NA', MAC:'NA', NS:'NA', CR:'NA', 
        POP_EH:'NA', POP_HWE:'NA', POP_MAF:'NA', POP_MAC:'NA', POP_NS:'NA', POP_CR:'NA', 
        GQ:'NA', gtDPmin:'NA', gtDPmax:'NA',
        MIN_SAMPLES_PER_POP:'NA', N_POPS_FAILING:'NA', PERC_POPS_FAILING:'NA',
        SAMPLE_MAX_MISSING:'NA',
    ]

    // Function to create cannonical filter string (reused later for population filters)
    def canon = { Map m, Map defaults = [:] ->
            def merged = defaults + m  // m overrides defaults
            merged.collect { k,v -> "${k}=${v == null ? 'NA' : v}" }
                    .sort()
                    .join(";")
        }

    // Function to pull filters from parameters by type
    def FiltersForType = { String vt ->
        def prefix = (vt=='snp'?'snp': vt=='indel'?'indel': vt=='invariant'?'inv': null)
        def p = { String k -> params.containsKey(k) ? params[k] : null }
        [
            QUAL_THR                     : p("${prefix}_global_qual"),
            DPmin                        : p("${prefix}_global_dp_min"),
            DIST_INDEL                   : p("${prefix}_global_dist_indel"),
            EH                           : p("${prefix}_global_eh"),
            HWE                          : p("${prefix}_global_hwe"),
            MAF                          : p("${prefix}_global_maf"),
            MAC                          : p("${prefix}_global_mac"),
            NS                           : p("${prefix}_global_min_samples"),
            CR                           : p("${prefix}_global_min_callrate"),
            POP_EH                       : p("${prefix}_pop_eh"),
            POP_HWE                      : p("${prefix}_pop_hwe"),
            POP_MAF                      : p("${prefix}_pop_maf"),
            POP_MAC                      : p("${prefix}_pop_mac"),
            POP_NS                       : p("${prefix}_pop_min_samples"),
            POP_CR                       : p("${prefix}_pop_min_callrate"),
            GQ                           : p("gt_qual"),
            gtDPmin                      : p("gt_dp_min"),
            gtDPmax                      : p("gt_dp_max"),
            MIN_SAMPLES_PER_POP          : p("min_samples_per_pop"),
            N_POPS_FAILING               : p("n_pops_failing"),
            PERC_POPS_FAILING            : p("perc_pops_failing"),
            SAMPLE_MAX_MISSING           : p("sample_max_missing"),
        ]
    }

    // Add filters in key:value format as an extra element
    ch_vcf_types
        .combine(ch_dp_bounds, by: 0)
        .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi, dpLo, dpHi ->
            def gf = FiltersForType(vt)
            gf.DPlower = dpLo
            gf.DPupper = dpHi
            tuple(vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi,
                    canon(gf, FILTER_DEFAULTS))
        }
        .set { ch_vcf_types_filters }

    // Create popmap tsv file for population based filtering
    ch_sample_pop
        .map { sample, pop -> "${sample}\t${pop}\n" }
        .collectFile(name: 'popmap.tsv', newLine: false)
        .set { ch_popmap }

    // Global sites filters 
    // This has 3 outputs: 1 = filtered VCF, = sitelist vcf tagged with QC (.out.tagged_sitelist), 3 = filtered sitelist vcf (.out.filtered_sitelist)
    FILTER_VCF (
        ch_vcf_types_filters,
	    ch_mask_bed_vcf,
        ch_popmap.first(),
        MERGE_CHUNK_MISSING.out.missing_summary
    )

    // Use counts file to remove those chunks which contain no variants
    FILTER_VCF.out.vcf
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi) }
        .set { ch_filtered_vcf }

     /*
        Create site filtering QC plots
    */
    FILTER_VCF.out.tagged_sitelist
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n ->
            tuple(variant_type, interval_hash, vcf, tbi)
        }
    .set { ch_vcfs_for_qc }

    // Create site histograms - uses the tags from the soft filtered vcf
    CREATE_FILTER_HIST(
        ch_vcfs_for_qc        
    )

    // QC plots for site histograms
    PLOT_VCF_FILTERS (
        CREATE_FILTER_HIST.out.hist.collect(),
        CREATE_FILTER_HIST.out.summary.collect(),
        "site_filters"
    )

     /*
        Create sitelist files for re-genotyping
    */
    
    FILTER_VCF.out.sitelist
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, counts_file ->
            def n = counts_file.text.trim() as Integer
            tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n)
        }
        .filter { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> n > 0 }
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, n -> tuple(variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi) }
        .set { ch_filtered_sites }

    // Merge sitelists back together by interval
    ch_filtered_sites.map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi-> tuple(interval_hash, vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_sitelists_to_merge }

    // Group all filtered sitelists by variant type and merge
    MERGE_FILTERED_SITELISTS (
        ch_sitelists_to_merge
    )

    // Reattach interval metadata
    MERGE_FILTERED_SITELISTS.out
        .join(ch_filtered_sites
        .map { variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi ->
            tuple(interval_hash, variant_type, interval_bed, bed_tbi)
        })
        .map { interval_hash, merged_vcf, merged_tbi, variant_type, interval_bed, bed_tbi ->
            tuple(variant_type, interval_hash, interval_bed, bed_tbi, merged_vcf, merged_tbi)
        }
        .set { ch_filtered_sites_merged }

    FILTER_VCF.out.samples_to_keep.first()
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }
   
    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_vcf = ch_filtered_vcf
    filtered_sitelist = ch_filtered_sites_merged
    sample_names_filt = ch_sample_names_filt
}