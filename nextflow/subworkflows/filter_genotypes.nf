/*
    Filter genotypes in .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_GT                                } from '../modules/filter_vcf_gt'
//include { PLOT_GT_FILTERS                              } from '../modules/plot_variant_qc'


workflow FILTER_GENOTYPES {

    take:
    ch_vcf
    ch_genome_indexed
    ch_mask_bed_vcf
    ch_sample_names

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
    FILTER_GENOTYPES (
        ch_vcf,
        ch_geno_filters
    )

    // plot genotype qc
    //PLOT_VARIANT_QC (
    //    FILTER_GENOTYPES.out.tables,
    //    ch_geno_filters
    //)
        
    emit: 
    vcf = FILTER_GENOTYPES.out.vcf

}