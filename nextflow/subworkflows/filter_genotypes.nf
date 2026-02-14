/*
    Filter .vcf files from GATK
*/

//// import modules
include { FILTER_VCF_GENOTYPES                         } from '../modules/filter_vcf_genotypes'
include { FILTER_VCF_MISSING                           } from '../modules/filter_vcf_missing'
include { PLOT_VCF_FILTERS as PLOT_GENOTYPE_FILTERS    } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'
include { CALC_CHUNK_MISSING                           } from '../modules/calc_chunk_missing'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing'

workflow FILTER_GENOTYPES {

    take:
    ch_genotyped_all

    main: 


    // Filter genotypes for quality - Set to missing genotype but retain GL/PL for probabilistic analyses
    FILTER_VCF_GENOTYPES (
        ch_genotyped_all
    )

    // TODO: FIlter_VCF_Genotypes outputs need to be named by chunk

    // QC plots for genotypes
    PLOT_GENOTYPE_FILTERS (
        FILTER_VCF_GENOTYPES.out.hist.collect(),
        FILTER_VCF_GENOTYPES.out.summary.collect(),
        "genotype_filters"
    )

    // Calculate per-sample missing data from each chunk
    CALC_CHUNK_MISSING (
        FILTER_VCF_GENOTYPES.out.vcf
    )

    // Merge per-sample missing data and DP histogram from all chunks into a single table
    MERGE_CHUNK_MISSING (
         CALC_CHUNK_MISSING.out.chunk_missing.map { variant_type, interval_hash, interval_bed, bed_tbi, missing -> missing }.collect()
    )

    // Filter for missing data
    FILTER_VCF_MISSING (
        FILTER_VCF_GENOTYPES.out.vcf,
        MERGE_CHUNK_MISSING.out.missing_summary.first()
    )

    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS (
        MERGE_CHUNK_MISSING.out.missing_summary,
        params.sample_max_missing
    )

    FILTER_VCF_MISSING.out.samples_to_keep
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }


    // Subset the merged vcf channels to each variant type for emission
    emit:
    vcf = FILTER_VCF_MISSING.out.vcf
    sample_names_filt = ch_sample_names_filt

}