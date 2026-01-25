/*
    Filter .vcf files from GATK
*/

//// import modules
include { CALC_CHUNK_MISSING                           } from '../modules/calc_chunk_missing'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing'
include { FILTER_VCF_GENOTYPES                         } from '../modules/filter_vcf_sites'
include { MERGE_VCFS as MERGE_FILTERED_VCFS            } from '../modules/merge_vcfs'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_GENOTYPES {

    take:
    ch_genotyped_all

    main: 


    // Filter genotypes for quality - Set to missing genotype but retain GL/PL for probabilistic analyses
    // Apply final missing data
    FILTER_VCF_GENOTYPES (
        ch_genotyped_all
    )


    // QC plots for sites and genotypes
    PLOT_VCF_FILTERS (
        FILTER_VCF_SITES.out.tables.collect()
    )

    // QC plots for sample missing data
    //PLOT_SAMPLE_FILTERS (
    //    MERGE_CHUNK_MISSING.out.missing_summary,
    //    params.sample_max_missing
    //)

    //FILTER_VCF.out.samples_to_keep
    //    .splitText( by: 1 )
    //    .unique()
    //    .set { ch_sample_names_filt }

        

    // Extract filtered sites only
    EXTRACT_VCF_SITES (
        MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='combined' }
    )   


    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_combined = ch_combined_filtered
    filtered_snps = ch_snp_filtered
    filtered_indels = ch_indel_filtered
    reports = VCF_STATS.out.vcfstats

}