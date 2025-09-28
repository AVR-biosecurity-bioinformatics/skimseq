/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_VCF                                   } from '../modules/filter_vcf'
include { MERGE_VCFS                                   } from '../modules/merge_vcfs'
include { MERGE_VCFS as MERGE_ALL                      } from '../modules/merge_vcfs'
include { VCF_STATS                                    } from '../modules/vcf_stats'
include { PLOT_SITE_FILTERS                            } from '../modules/plot_site_filters'
include { PLOT_GT_FILTERS                              } from '../modules/plot_gt_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_SITES {

    take:
    ch_vcf
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

    // Create joint channel for next step
    Channel
    .from(SNP_FILTERS, INDEL_FILTERS, INV_FILTERS)
    .map { f -> tuple(f.type as String, f) }
    .set { ch_type_filters }     // emits: val(type), val(filters)

    // For each input VCF, run FILTER_VCF three times (snp/indel/invariant)
    ch_vcf
        .cross( ch_type_filters ) 
        .map { vcf, vcf_tbi, type, filters -> tuple(vcf, vcf_tbi, type, filters, ch_mask_bed_vcf) }
        | FILTER_VCF

    // Group all filtered VCFs by variant type and merge
    FILTER_VCF.out.vcf
    .map { type, vcf, tbi -> tuple(type, vcf) }
    .groupTuple()
    | MERGE_VCFS

    // plot variant qc
    // NOT WORKING WITH NEW MAP APPROACH
   //PLOT_SITE_FILTERS (
   //     FILTER_SNPS.out.tables,
   //     FILTER_INDELS.out.tables,
   //     FILTER_INVARIANT.out.tables,
   //     ch_snp_filters,
   //     ch_indel_filters,
   //     ch_inv_filters
   // )

    // Gather all per-type merged VCFs, then merge all into one
    MERGE_VCFS.out.vcf
        .map { type, vcf, tbi -> vcf }
        .collect() 
        .map { vcf_list -> tuple('all', vcf_list) }
        | MERGE_ALL
        .set { ch_merged_all }

    // Calculate VCF statistics
    VCF_STATS (
         ch_merged_all.map { sample, vcf, tbi -> [ vcf, tbi ] },
         ch_genome_indexed,
         ch_sample_names
    )
        
    emit: 
    filtered_all = ch_merged_all.map { type, vcf, tbi -> [ vcf, tbi ] }
    filtered_snps = MERGE_VCFS.out.vcf.filter{ it[0]=='snp' }.map{ _, vcf, tbi -> [vcf,tbi] }
    filtered_indels = MERGE_VCFS.out.vcf.filter{ it[0]=='indel' }.map{ _, vcf, tbi -> [vcf,tbi] }
    reports = VCF_STATS.out.vcfstats

}