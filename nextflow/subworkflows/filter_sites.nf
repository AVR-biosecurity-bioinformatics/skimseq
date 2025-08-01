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


workflow FILTER_SITES {

    take:
    ch_vcf
    ch_genome_indexed
    ch_mask_bed_vcf
    ch_sample_names

    main: 

    // collect SNP filtering parameters into a single list
    Channel.of(
        "snp",
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
        params.snp_max_missing,
        params.snp_custom_flags
    )
    .collect( sort: false )
    .set { ch_snp_filters }

    // collect indel filtering parameters into a single list
    Channel.of (
        "indel",
        params.indel_qd,          
        params.indel_qual,        
        "NA",         
        params.indel_fs,          
        "NA",       
        "NA",       
        params.indel_rprs,        
        params.indel_maf,      
        params.indel_mac,            
        params.indel_eh,          
        params.indel_dp_min,      
        params.indel_dp_max,
        params.indel_max_missing,
        params.indel_custom_flags
    )
    .collect ( sort: false )
    .set { ch_indel_filters }

    // collect invariant filtering parameters into a single list
    Channel.of (     
        "invariant",
        "NA",          
        "NA",       
        "NA",         
        "NA",          
        "NA",       
        "NA",       
        "NA",    
        "NA",     
        "NA",            
        "NA",          
        params.inv_dp_min,      
        params.inv_dp_max,
        params.inv_max_missing,
        params.inv_custom_flags,
    )
    .collect ( sort: false )
    .set { ch_inv_filters }

    // filter SNPs
    FILTER_SNPS (
        ch_vcf,
        ch_snp_filters,
        ch_mask_bed_vcf
    )

    // filter indels
    FILTER_INDELS (
        ch_vcf,
        ch_indel_filters,
        ch_mask_bed_vcf
    )

    // filter indels
    FILTER_INVARIANT (
        ch_vcf,
        ch_inv_filters,
        ch_mask_bed_vcf
    )

    // plot variant qc
    PLOT_SITE_FILTERS (
        FILTER_SNPS.out.tables,
        FILTER_INDELS.out.tables,
        FILTER_INVARIANT.out.tables,
        ch_snp_filters,
        ch_indel_filters,
        ch_inv_filters
    )

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
         MERGE_FILTERED.out.vcf,
         ch_genome_indexed,
         ch_sample_names
    )

    // Create reports channel for multiqc
    VCF_STATS.out.vcfstats
        .set { ch_reports}
        
    emit: 
    filtered_merged = MERGE_FILTERED.out.vcf
    filtered_snps = FILTER_SNPS.out.vcf
    filtered_indels = FILTER_INDELS.out.vcf
    reports = ch_reports

}