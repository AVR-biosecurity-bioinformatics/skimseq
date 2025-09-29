/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_VCF                                   } from '../modules/filter_vcf'
include { MERGE_VCFS as MERGE_FILTERED_VCFS            } from '../modules/merge_vcfs'
include { VCF_STATS                                    } from '../modules/vcf_stats'
include { PLOT_SITE_FILTERS                            } from '../modules/plot_site_filters'
include { PLOT_GT_FILTERS                              } from '../modules/plot_gt_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_SITES {

    take:
    ch_vcfs
    ch_genome_indexed
    ch_mask_bed_vcf
    ch_sample_names

    main: 

    // collect SNP filtering parameters into a map
    // TODO: how to handle custom flags?
    def SNP_FILTERS = [
        type: 'snp',
        qd: params.snp_qd,
        qual: params.snp_qual,
        sor: params.snp_sor,
        fs: params.snp_fs,
        mq: params.snp_mq,
        mqrs: params.snp_mqrs,
        rprs: params.snp_rprs,
        maf: params.snp_maf,
        mac: params.snp_mac,
        eh: params.snp_eh,
        dp_min: params.snp_dp_min,
        dp_max: params.snp_dp_max,
        max_missing: params.snp_max_missing,
        gq: params.gt_qual,
        gt_dp_min: params.gt_dp_min,
        gt_dp_max: params.gt_dp_max
    ]

    // collect indel filtering parameters into a map
    def INDEL_FILTERS = [
        type: 'indel',
        qd: params.indel_qd,
        qual: params.indel_qual,
        fs: params.indel_fs,
        rprs: params.indel_rprs,
        maf: params.indel_maf,
        mac: params.indel_mac,
        eh: params.indel_eh,
        dp_min: params.indel_dp_min,
        dp_max: params.indel_dp_max,
        max_missing: params.indel_max_missing,
        gq: params.gt_qual,
        gt_dp_min: params.gt_dp_min,
        gt_dp_max: params.gt_dp_max
    ]

    // collect invariant filtering parameters into a single list
    def INV_FILTERS = [
        type: 'invariant',
        dp_min: params.inv_dp_min,
        dp_max: params.inv_dp_max,
        max_missing: params.inv_max_missing,
        gq: params.gt_qual,
        gt_dp_min: params.gt_dp_min,
        gt_dp_max: params.gt_dp_max
    ]

    // Create joint of all 3 variant type filters
    Channel
        .from(SNP_FILTERS, INDEL_FILTERS, INV_FILTERS)
        .map { f -> tuple(f.type as String, f) }
        .set { ch_type_filters }
    
    // For each input VCF, combine with type_filters and run FILTER_VCF for each type
    FILTER_VCF (
        ch_vcfs.combine( ch_type_filters ),
	    ch_mask_bed_vcf
    )

    // Create a channel of all 3 variant types + all together for merging
    FILTER_VCF.out.vcf
        .concat(FILTER_VCF.out.vcf.map { type, vcf, tbi -> tuple('combined', vcf, tbi) })
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    // Group all filtered VCFs by variant type and merge
    MERGE_FILTERED_VCFS (
        ch_vcf_to_merge
    )
      
    // Extract merged variant type vcfs into convenient channels
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='combined' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_combined_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='snp' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_snp_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='indel' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_indel_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='invariant' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_invariant_filtered }

    // Variant QC plots
    // NOT WORKING WITH NEW MAP APPROACH

    // plot site qc
    PLOT_SITE_FILTERS (
        FILTER_VCF.out.tables
    )

    // plot genotype qc
    //PLOT_GT_FILTERS (
    //    FILTER_VCF_GT.out.tables,
    //    ch_geno_filters
    //)

    // plot samples qc
    //PLOT_SAMPLE_FILTERS (
    //    FILTER_VCF_SAMPLES.out.tables,
    //    params.sample_max_missing
    //)

    // Calculate VCF statistics
   VCF_STATS (
        ch_combined_filtered,
        ch_genome_indexed,
        ch_sample_names
    )
        
    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_combined = ch_combined_filtered
    filtered_snps = ch_snp_filtered
    filtered_indels = ch_indel_filtered
    reports = VCF_STATS.out.vcfstats

}