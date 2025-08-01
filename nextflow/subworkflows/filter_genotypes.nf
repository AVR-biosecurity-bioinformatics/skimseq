/*
    Filter genotypes in .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_GT                                } from '../modules/filter_vcf_gt'
include { PLOT_GT_FILTERS                              } from '../modules/plot_gt_filters'


workflow FILTER_GENOTYPES {

    take:
    ch_vcf

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
        
    emit: 
    vcf = FILTER_VCF_GT.out.vcf

}