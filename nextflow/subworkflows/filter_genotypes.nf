/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_GENOTYPES                         } from '../modules/filter_vcf_genotypes'
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
   // PLOT_VCF_FILTERS (
    //    FILTER_VCF_SITES.out.tables.collect()
    //)

    // QC plots for sample missing data
    //PLOT_SAMPLE_FILTERS (
    //    MERGE_CHUNK_MISSING.out.missing_summary,
    //    params.sample_max_missing
    //)

    //FILTER_VCF.out.samples_to_keep
    //    .splitText( by: 1 )
    //    .unique()
    //    .set { ch_sample_names_filt }



    // Subset the merged vcf channels to each variant type for emission
    emit:
    vcf = FILTER_VCF_GENOTYPES.out.vcf

}