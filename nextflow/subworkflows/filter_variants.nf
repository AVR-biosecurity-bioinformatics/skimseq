/*
    Filter .vcf files from GATK
*/

//// import modules
include { CALC_CHUNK_DP                                } from '../modules/calc_chunk_dp'
include { MERGE_CHUNK_DP                               } from '../modules/merge_chunk_dp'
include { FILTER_VCF_SITES as FILTER_VCF_SITES_GLOBAL  } from '../modules/filter_vcf_sites'
include { FILTER_VCF_SITES as FILTER_VCF_SITES_POP     } from '../modules/filter_vcf_sites'
include { INTERSECT_FILTERED_SITES                     } from '../modules/intersect_filtered_sites'
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

    /*
        Function definitions
    */

    // Set up default emtpy filters
    def FILTER_DEFAULTS = [
        QUAL_THR:'NA', DPlower:'NA', PCT_LOW:'NA', DPupper:'NA', DIST_INDEL:'NA',
        EH:'NA', HWE:'NA', MAF:'NA', MAC:'NA', NS:'NA', CR:'NA'
    ]

   // Function to create cannonical filter string
   def canon = { Map m, Map defaults = [:] ->
        def merged = defaults + m  // m overrides defaults
        merged.collect { k,v -> "${k}=${v == null ? 'NA' : v}" }
                .sort()
                .join(";")
    }

    // Duplicate vcf files by variant type
    def variant_types = ['snp', 'indel']
        if( params.output_invariant ) {
        variant_types << 'invariant'
    }
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

    /*
        Global filters for sites
    */

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
  ch_sample_pop
    .map { sample, pop -> tuple(pop, sample) }   // key by pop
    .groupTuple()                                // -> (pop, [sample1, sample2, ...])
    .map { pop, samples -> tuple(pop, [samples: samples.unique().sort()]) }
    .set { ch_sample_list_pop }

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
        .map { vt, ih, interval_bed, bed_tbi, global_vcf, global_tbi, counts ->
        tuple([vt, ih], interval_bed, bed_tbi, global_vcf, global_tbi)
        }
        .join(  FILTER_VCF_SITES_POP.out.vcf
            .map { vt, pop, ih, interval_bed, bed_tbi, pop_vcf, pop_tbi, counts ->
            tuple([vt, ih], pop_vcf)
            }
            .groupTuple()
        )    
    .map { key, interval_bed, bed_tbi, global_vcf, global_tbi, pop_vcfs ->
      def (vt, ih) = key
      tuple(vt, ih, interval_bed, bed_tbi, global_vcf, pop_vcfs)
    }
    .set { ch_sitelists_to_intersect }

    // intersect global filter with per-pop, keeping only those failing in >n populations
    // This function makes the QC histograms too
    INTERSECT_FILTERED_SITES (
        ch_sitelists_to_intersect
    )

    // Create site histograms - uses the tages from the soft filtered vcf
    //CREATE_FILTER_HISTS(
    //    
    //)

    // QC plots for site histograms
    //PLOT_VCF_FILTERS (
    //    FILTER_VCF_SITES_GLOBAL.out.hist.collect(),
    //    FILTER_VCF_SITES_GLOBAL.out.summary.collect(),
    //    "site_filters"
    //)


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