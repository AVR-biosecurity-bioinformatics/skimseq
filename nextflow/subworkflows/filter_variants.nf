/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_INDELS                                                 } from '../modules/filter_indels'
include { FILTER_SNPS                                                 } from '../modules/filter_snps'
include { MERGE_FILTERED                                                 } from '../modules/merge_filtered'



workflow FILTER_VARIANTS {

    take:
    ch_vcf

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


    FILTER_SNPS (
        ch_vcf,
        ch_snp_filters,
        params.max_missing
    )
    
    FILTER_INDELS (
        ch_vcf
    )

    MERGE_FILTERED (
        FILTER_SNPS.out.vcf,
        FILTER_INDELS.out.vcf
    )



    emit: 
    FILTER_SNPS.out.vcf


}