/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_INDELS                                            } from '../modules/filter_indels'
include { FILTER_SNPS                                              } from '../modules/filter_snps'
include { FILTER_INVARIANT                                         } from '../modules/filter_invariant'
include { MERGE_FILTERED                                           } from '../modules/merge_filtered'
include { VCF_STATS                                                } from '../modules/vcf_stats'
include { PLOT_VARIANT_QC                                          } from '../modules/plot_variant_qc'



workflow FILTER_VARIANTS {

    take:
    ch_vcf
    ch_genome_indexed
    ch_mask_bed_vcf

    main: 

    // collect generic genotype filtering parameters into a single list
    // These are used for all variant types
    Channel.of(
        params.max_nocall,
        params.max_missing,
        params.gt_qual,        
        params.gt_dp_min,         
        params.gt_dp_max
    )
    .collect( sort: false )
    .set { ch_geno_filters }

    // collect SNP filtering parameters into a single list
    Channel.of(
        params.snp_qd,          
        params.snp_qual,        
        params.snp_sor,         
        params.snp_fs,          
        params.snp_mq,          
        params.snp_mqrs,        
        params.snp_rprs,        
        params.snp_maf,      
        params.snp_mac,            
        params.snp_eh,          
        params.snp_dp_min,      
        params.snp_dp_max,      
        params.snp_custom_flags,
    )
    .collect( sort: false )
    .set { ch_snp_filters }

    // filter SNPs
    FILTER_SNPS (
        ch_vcf,
        ch_snp_filters,
        ch_geno_filters,
        ch_mask_bed_vcf
    )
    
    // collect indel filtering parameters into a single list
    Channel.of (
        params.indel_qd,          
        params.indel_qual,        
        params.indel_fs,          
        params.indel_rprs,        
        params.indel_maf,      
        params.indel_mac,            
        params.indel_eh,          
        params.indel_dp_min,      
        params.indel_dp_max,      
        params.indel_custom_flags,
    )
    .collect ( sort: false )
    .set { ch_indel_filters }

    // filter indels
    FILTER_INDELS (
        ch_vcf,
        ch_indel_filters,
        ch_geno_filters,
        ch_mask_bed_vcf
    )

    // collect invariant filtering parameters into a single list
    Channel.of (     
        params.inv_dp_min,      
        params.inv_dp_max,      
        params.inv_custom_flags,
    )
    .collect ( sort: false )
    .set { ch_inv_filters }

    // filter indels
    FILTER_INVARIANT (
        ch_vcf,
        ch_inv_filters,
        ch_geno_filters,
        ch_mask_bed_vcf
    )

    // plot variant qc
    PLOT_VARIANT_QC (
        FILTER_SNPS.out.tables,
        FILTER_INDELS.out.tables,
        FILTER_INVARIANT.out.tables,
        ch_snp_filters,
        ch_indel_filters,
        ch_inv_filters,
        params.max_missing
    )

    // merge filtered SNPs and indels together into one file
    MERGE_FILTERED (
        FILTER_SNPS.out.vcf,
        FILTER_INDELS.out.vcf,
        FILTER_INVARIANT.out.vcf
    )

    // Calculate VCF statistics
    VCF_STATS (
         MERGE_FILTERED.out.vcf,
         ch_genome_indexed
    )

    // Create reports channel for multiqc
    VCF_STATS.out.vcfstats
        .set { ch_reports}
        
    emit: 
    filtered_vcf = MERGE_FILTERED.out.vcf
    reports = ch_reports

}