/*
    Filter .vcf files from GATK
*/

//// import modules
include { EXTRACT_VCF_SITES                            } from '../modules/extract_vcf_sites'
include { COUNT_VCF_RECORDS                            } from '../modules/count_vcf_records'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_FILT    } from '../modules/create_interval_chunks'
include { SUBSET_VCF_TO_SITES                          } from '../modules/subset_vcf_to_sites'
include { CALC_DP_BOUNDS                               } from '../modules/calc_dp_bounds'
include { FILTER_VCF_SITES as FILTER_VCF_SITES_GLOBAL  } from '../modules/filter_vcf_sites'
include { FILTER_VCF_SITES as FILTER_VCF_SITES_POP     } from '../modules/filter_vcf_sites'
include { INTERSECT_FILTERED_SITES                     } from '../modules/intersect_filtered_sites'
include { MERGE_VCFS as MERGE_UNFILTERED_SITELISTS     } from '../modules/merge_vcfs'
include { MERGE_VCFS as MERGE_FILTERED_SITELISTS       } from '../modules/merge_vcfs'
include { CREATE_FILTER_HIST                           } from '../modules/create_filter_hist'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'

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

    /*
        Create new interval chunks based on the number of called variants
    */

    // Extract a sites-only VCF (drop genotypes) of unfiltered sites for faster computation
    EXTRACT_VCF_SITES(
        ch_vcfs.map { interval_hash, interval_bed, bed_tbi, vcf, tbi -> tuple(interval_hash, vcf, tbi) }
    )

    // Merge per-chunk unfiltered sites lists into a single vcf
    EXTRACT_VCF_SITES.out.vcf
        .map { interval_hash, vcf, tbi -> tuple('unfiltered_sitelist', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_unfiltered_sitelist_to_merge }

    MERGE_UNFILTERED_SITELISTS (
        ch_unfiltered_sitelist_to_merge
    )

    // Count number of VCF records for chunking, this also outputs a depth histogram used later for calculating DP bounds
    COUNT_VCF_RECORDS (
        MERGE_UNFILTERED_SITELISTS.out.vcf,
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_mask_bed_vcf
    )

    // Create new interval chunks from the merged unfiltered sitelist
    CREATE_INTERVAL_CHUNKS_FILT (
        COUNT_VCF_RECORDS.out.counts,
        ch_genome_indexed,
        ch_include_bed.first(),
        params.filt_var_per_chunk,
        0,
        "false",
        "false"
    )

    // Parse interval outputs
    CREATE_INTERVAL_CHUNKS_FILT.out.interval_bed
            .flatMap { sample, beds, tbis  ->
                // normalize to a list for cases where there are only 1 bed output for a sample
                def bedList = (beds instanceof List) ? beds : [beds]
                def tbiList = (tbis instanceof List) ? tbis : [tbis]

                assert bedList.size() == tbiList.size() :
                "Mismatch for ${sample}: beds=${bedList.size()} tbis=${tbiList.size()}"

                // emit one tuple per bed file
                (0..<bedList.size()).collect { i ->
                    def bed = bedList[i] as Path
                    def bed_tbi = tbiList[i]
                    def base = bed.getFileName().toString()
                    base = base.replaceFirst(/\.gz$/, '')
                    base = base.replaceFirst(/\.bed$/, '')
                    def interval_hash = base.startsWith('_') ? base.substring(1) : base
                    tuple(interval_hash, bed, bed_tbi)
                }
            }
            .set { ch_interval_bed_filt }

    // Create a file with both the original genotypes and the interval bed
    ch_vcfs
        .map { interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi -> tuple(interval_hash, vcf, vcf_tbi) }
        .toSortedList { a, b -> a[0] <=> b[0] }          // sort by interval_hash to ensure its deterministicly ordered
        .map { rows ->
            def vcfs = rows.collect { it[1] }
            def tbis = rows.collect { it[2] }
            tuple(vcfs, tbis)
        }
        .set { ch_vcfs_list }

    ch_interval_bed_filt
        .combine(ch_vcfs_list)
        .map { interval_hash, bed, bed_tbi, vcf_list, vcf_tbi_list -> tuple("all", interval_hash, bed, bed_tbi, vcf_list, vcf_tbi_list ) }
        .set { ch_sitelist_with_all_vcfs }


    // Extract the sites in each new chunks from the original (with genotype) vcf file
    SUBSET_VCF_TO_SITES(
        ch_sitelist_with_all_vcfs
    )

    SUBSET_VCF_TO_SITES.out.vcf
        .map { variant_type, interval_hash, bed, bed_tbi, vcf, vcf_tbi -> tuple(interval_hash, bed, bed_tbi, vcf, vcf_tbi  ) }
        .set { ch_vcfs_rechunked }

    /*
        Calculate depth filters
    */

    // Calculate chunk_dp by varaint type using the dphist from COUNT_VCF_RECORDS
    Channel.of(*variant_types)
        .combine(COUNT_VCF_RECORDS.out.dphist)
        .map { variant_type, name, dphist -> tuple(variant_type, dphist ) }
        .groupTuple()
        .map { vt, files ->
            def lo
            def hi
            switch(vt) {
                case 'snp':
                lo = params.snp_global_dp_lower_perc
                hi = params.snp_global_dp_upper_perc
                break
                case 'indel':
                lo = params.indel_global_dp_lower_perc
                hi = params.indel_global_dp_upper_perc
                break
                case 'invariant':
                lo = params.inv_global_dp_lower_perc
                hi = params.inv_global_dp_upper_perc
                break
                default:
                lo = null; hi = null
            }
            tuple(vt, files, lo, hi)
        }
        .set { ch_chunk_dp_grouped }

    CALC_DP_BOUNDS (
        ch_chunk_dp_grouped
    )

    CALC_DP_BOUNDS.out.dp_bounds.map { vt, f ->
        def lines = f.readLines()
        def hdr = lines[0].split('\t')
        def row = lines[1].split('\t')
        def m = [hdr, row].transpose().collectEntries { k,v -> [(k): v] }
        tuple(vt, m.DPlower as Integer, m.DPupper as Integer)
    }
    .set { ch_dp_bounds }

    /*
        Global filters for sites
    */

    // Set up default emtpy filters
        def FILTER_DEFAULTS = [
            QUAL_THR:'NA', DPlower:'NA', PCT_LOW:'NA', DPupper:'NA', DIST_INDEL:'NA',
            EH:'NA', HWE:'NA', MAF:'NA', MAC:'NA', NS:'NA', CR:'NA'
        ]

    // Function to create cannonical filter string (reused later for population filters)
    def canon = { Map m, Map defaults = [:] ->
            def merged = defaults + m  // m overrides defaults
            merged.collect { k,v -> "${k}=${v == null ? 'NA' : v}" }
                    .sort()
                    .join(";")
        }
    // Function to pull global filters from parameters by type
    def GlobalFiltersForType = { String vt ->
        def prefix = (vt=='snp'?'snp': vt=='indel'?'indel': vt=='invariant'?'inv': null)
        def p = { String k -> params.containsKey(k) ? params[k] : null }
        [
            QUAL_THR  : p("${prefix}_global_qual"),
            DPmin     : p("${prefix}_global_dp_min"),
            DIST_INDEL: p("${prefix}_global_dist_indel"),
        ]
    }

    // Duplicate vcf files by variant type
    Channel.of(*variant_types)
        .combine(ch_vcfs_rechunked)
        .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi ->
        tuple(vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi)
        }
    .set { ch_vcf_types }

    ch_sample_names
        .collect()
        .map { it.unique().sort() } 
        .map { [samples: it] }
        .set{ ch_sample_list_global }

    // Create channel of VCFs by variant type, with filter string as an extra element
    ch_vcf_types
        .join(ch_dp_bounds)
        .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi, dpLo, dpHi ->
            def gf = GlobalFiltersForType(vt)
            gf.DPlower = dpLo
            gf.DPupper = dpHi
            tuple(vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi,
                    canon(gf, FILTER_DEFAULTS))
        }
        .combine(ch_sample_list_global)   // attaches samples file to every tuple
        .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi, filter_kv, sample_names ->
            tuple(vt, "global", interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi, filter_kv, sample_names )
        }
        .set { ch_vcf_types_global }

    // Global sites filters 
    FILTER_VCF_SITES_GLOBAL (
        ch_vcf_types_global,
	    ch_mask_bed_vcf
    )

    /*
        Per-population filters for sites
    */

    // Function to pull per-population filters from parameters by type
    def PopFiltersForType = { String vt ->
        def prefix = (vt=='snp'?'snp': vt=='indel'?'indel': vt=='invariant'?'inv': null)
        def p = { String k -> params.containsKey(k) ? params[k] : null }
        [
            EH        : p("${prefix}_pop_eh"),
            HWE       : p("${prefix}_pop_hwe"),
            MAF       : p("${prefix}_pop_maf"),
            MAC       : p("${prefix}_pop_mac"),
            NS        : p("${prefix}_pop_min_samples"),
            CR        : p("${prefix}_pop_min_callrate"),
        ]
    }

  // Create seperate channels that have all sample names grouped by pop
  // Exclude any below min samples per pop
  ch_sample_pop
    .map { sample, pop -> tuple(pop, sample) }
    .groupTuple()
    .map { pop, samples -> tuple(pop, samples.unique().sort()) }
    .filter { pop, samples -> samples.size() >= (params.min_samples_per_pop as Integer) }
    .map { pop, samples -> tuple(pop, [samples: samples]) }
    .set {ch_sample_list_pop }

    // Create channel of VCFs by variant type, with filter string as an extra element
    ch_vcf_types
    .combine(ch_sample_list_pop)   // (vt, ih, bed..., vcf..., pop, samples)
    .map { vt, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi, pop, samples ->
        tuple(vt, pop, interval_hash, interval_bed, bed_tbi, vcf, vcf_tbi,
            canon(PopFiltersForType(vt), FILTER_DEFAULTS),
            samples)
    }
    .set { ch_vcf_types_pop }

    // Per-population site filters
    FILTER_VCF_SITES_POP (
        ch_vcf_types_pop,
	    ch_mask_bed_vcf
    )

    /*
        Intersect global and per-population filters
    */

    FILTER_VCF_SITES_GLOBAL.out.vcf
    .map { vt, ih, interval_bed, bed_tbi, global_vcf, global_tbi ->
      tuple([vt, ih], interval_bed, bed_tbi, global_vcf, global_tbi)
    }
    .set{ ch_global_keyed }

  FILTER_VCF_SITES_POP.out.vcf
    .map { vt, ih, interval_bed, bed_tbi, pop_vcf, pop_tbi ->
      tuple([vt, ih], tuple(pop_vcf, pop_tbi))
    }
    .groupTuple()
    .set{ ch_pop_grouped }

  ch_global_keyed
    .join(ch_pop_grouped)
    .map { key, interval_bed, bed_tbi, global_vcf, global_tbi, pop_pairs ->
      def (vt, ih) = key
      def pop_vcfs = pop_pairs.collect{ it[0] }
      def pop_tbis = pop_pairs.collect{ it[1] }
      tuple(vt, ih, interval_bed, bed_tbi, global_vcf, global_tbi, pop_vcfs, pop_tbis)
    }
    .set { ch_sitelists_to_intersect }

    // intersect global filter with per-pop, keeping only those failing in >n populations
    // This function makes the QC histograms too
    INTERSECT_FILTERED_SITES (
        ch_sitelists_to_intersect,
        params.n_pops_failing,
        params.perc_pops_failing
    )

     /*
        Create site filtering QC plots
    */
    FILTER_VCF_SITES_GLOBAL.out.vcf
        .concat(FILTER_VCF_SITES_POP.out.vcf)    
        .map { vt, ih, interval_bed, bed_tbi, vcf, tbi ->
            tuple(vt, ih, vcf, tbi)
        }
    .set { ch_vcfs_for_qc }

    // Create site histograms - uses the tages from the soft filtered vcf
    CREATE_FILTER_HIST(
        ch_vcfs_for_qc        
    )

    // QC plots for site histograms
    PLOT_VCF_FILTERS (
        CREATE_FILTER_HIST.out.hist.collect(),
        CREATE_FILTER_HIST.out.summary.collect(),
        "site_filters"
    )

    // Use counts file to remove those chunks which contain no variants
    INTERSECT_FILTERED_SITES.out.vcf
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