/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_SITES as FILTER_SNPS                          } from '../modules/filter_vcf_sites'
include { FILTER_VCF_SITES as FILTER_INDELS                        } from '../modules/filter_vcf_sites'
include { FILTER_VCF_SITES as FILTER_INVARIANT                     } from '../modules/filter_vcf_sites'
include { MERGE_VCFS as MERGE_FILTERED                             } from '../modules/merge_vcfs'
include { VCF_STATS                                                } from '../modules/vcf_stats'
include { PLOT_SITE_FILTERS                                        } from '../modules/plot_site_filters'
include { FILTER_VCF_GT                                } from '../modules/filter_vcf_gt'
include { PLOT_GT_FILTERS                              } from '../modules/plot_gt_filters'
include { FILTER_VCF_SAMPLES                               } from '../modules/filter_vcf_samples'
include { PLOT_SAMPLE_FILTERS                              } from '../modules/plot_sample_filters'

workflow FILTER_SITES {

    take:
    ch_vcf
    ch_genome_indexed
    ch_mask_bed_vcf

    main: 

    // collect generic genotype filtering parameters into a single list
    // These are used for all variant types
    Channel.of(
        params.gt_qual,        
        params.gt_dp_min,         
        params.gt_dp_max
    )
    .collect( sort: false )
    .set { ch_geno_filters }

    // filter genotypes
    FILTER_VCF_GT (
        ch_vcf,
        ch_geno_filters
    )

    // plot genotype qc
    PLOT_GT_FILTERS (
        FILTER_VCF_GT.out.tables,
        ch_geno_filters
    )

    // filter samples
    FILTER_VCF_SAMPLES (
        FILTER_VCF_GT.out.vcf,
        params.sample_max_missing
    )

    FILTER_VCF_SAMPLES.out.samples_to_keep
        .splitText( by: 1 )
        .set { ch_sample_names }

    // plot samples qc
    PLOT_SAMPLE_FILTERS (
        FILTER_VCF_SAMPLES.out.tables,
        params.sample_max_missing
    )

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
        max_missing: params.snp_max_missing
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
        max_missing: params.indel_max_missing
    ]

    // collect invariant filtering parameters into a single list
    def INV_FILTERS = [
        type: 'invariant',
        dp_min: params.inv_dp_min,
        dp_max: params.inv_dp_max,
        max_missing: params.inv_max_missing
    ]

    // Create value channels
    Channel.value(SNP_FILTERS).set { ch_snp_filters }
    Channel.value(INDEL_FILTERS).set { ch_indel_filters }
    Channel.value(INV_FILTERS).set { ch_inv_filters }


    // filter SNPs
    FILTER_SNPS (
        FILTER_VCF_SAMPLES.out.vcf,
        "snp",
        ch_snp_filters,
        ch_mask_bed_vcf
    )

    // filter indels
    FILTER_INDELS (
        FILTER_VCF_SAMPLES.out.vcf,
        "indel",
        ch_indel_filters,
        ch_mask_bed_vcf
    )

    // filter indels
    FILTER_INVARIANT (
        FILTER_VCF_SAMPLES.out.vcf,
        "invariant",
        ch_inv_filters,
        ch_mask_bed_vcf
    )

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

    // Create channel of VCFs to merge
    FILTER_SNPS.out.vcf
        .mix(FILTER_INDELS.out.vcf, FILTER_INVARIANT.out.vcf)
        .collect(flat: false)
 	    .map { it.transpose() }
        .set { ch_vcfs }

    // merge filtered SNPs and indels together into one file
    // TODO: change the output name to the project name
    MERGE_FILTERED (
        ch_vcfs,
        "merged"
    )

    // Calculate VCF statistics
    VCF_STATS (
         MERGE_FILTERED.out.vcf.map { sample, vcf, tbi -> [ vcf, tbi ] },
         ch_genome_indexed,
         ch_sample_names
    )

    // Create reports channel for multiqc
    VCF_STATS.out.vcfstats
        .set { ch_reports}
        
    emit: 
    filtered_merged = MERGE_FILTERED.out.vcf.map { sample, vcf, tbi -> [ vcf, tbi ] }
    filtered_snps = FILTER_SNPS.out.vcf
    filtered_indels = FILTER_INDELS.out.vcf
    reports = ch_reports

}