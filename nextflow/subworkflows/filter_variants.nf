/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_INDELS                                                 } from '../modules/filter_indels'
include { FILTER_SNPS                                                 } from '../modules/filter_snps'
include { MERGE_FILTERED                                                 } from '../modules/merge_filtered'
include { VCF_STATS                                                } from '../modules/vcf_stats'



workflow FILTER_VARIANTS {

    take:
    ch_vcf
    ch_genome_indexed
    ch_mask_bed_vcf

    main: 

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
        params.max_missing,
        ch_mask_bed_vcf
    )
    
    // collect indel filtering parameters into a single list
    Channel.of (
        params.indel_qd,          
        params.indel_qual,        
        params.indel_fs,          
        params.indel_rprs,        
        params.indel_maf,         
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
        params.max_missing,
        ch_mask_bed_vcf
    )

    // merge filtered SNPs and indels together into one file
    MERGE_FILTERED (
        FILTER_SNPS.out.vcf,
        FILTER_INDELS.out.vcf
    )

    // Calculate VCF statistics
    VCF_STATS (
         MERGE_FILTERED.out.vcf,
         ch_genome_indexed
    )

    emit: 
    MERGE_FILTERED.out.vcf


}