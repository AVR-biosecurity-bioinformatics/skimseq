/*
    Filter samples in .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_SAMPLES                               } from '../modules/filter_vcf_samples'
//include { PLOT_SAMPLE_FILTERS                              } from '../modules/plot_sample_filters'


workflow FILTER_SAMPLES {

    take:
    ch_vcf

    main: 

    // collect generic genotype filtering parameters into a single list
    // These are used for all variant types
    Channel.of(
        params.sample_missing
    )
    .collect( sort: false )
    .set { ch_sample_filters }

    // filter samples
    FILTER_VCF_SAMPLES (
        ch_vcf,
        ch_sample_filters
    )

    // plot samples qc
    //PLOT_SAMPLE_FILTERS (
    //    FILTER_VCF_SAMPLES.out.tables,
    //    ch_sample_filters
    //)
        
    emit: 
    vcf = FILTER_VCF_SAMPLES.out.vcf

}