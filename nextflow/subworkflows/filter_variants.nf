/*
    Filter .vcf files from GATK
*/

//// import modules
include { CALC_DATASET_FILTERS                         } from '../modules/calc_dataset_filters'
include { FILTER_VCF                                   } from '../modules/filter_vcf'
include { MERGE_VCFS as MERGE_FILTERED_VCFS            } from '../modules/merge_vcfs'
include { VCF_STATS                                    } from '../modules/vcf_stats'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_VARIANTS {

    take:
    ch_vcfs
    ch_merged_vcf
    ch_genome_indexed
    ch_mask_bed_vcf

    main: 

    // Caclulate dataset-wide filters
    CALC_DATASET_FILTERS (
        ch_merged_vcf
    )

    // For each input VCF, combine with type to make 3 copies for each type, then run FILTER_VCF for each type
    FILTER_VCF (
        ch_vcfs.combine( channel.of("snp", "indel", "invariant") ),
	    ch_mask_bed_vcf,
        CALC_DATASET_FILTERS.out.missing_summary.first(),
        CALC_DATASET_FILTERS.out.dp_summary.first()
    )

    // Create a channel of all 3 variant types + all together for merging
    FILTER_VCF.out.vcf
        .concat(FILTER_VCF.out.vcf.map { type, vcf, tbi -> tuple('combined', vcf, tbi) })
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    // Group all filtered VCFs by variant type and merge
    MERGE_FILTERED_VCFS (
        ch_vcf_to_merge
    )
      
    // Extract merged variant type vcfs into convenient channels
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='combined' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_combined_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='snp' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_snp_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='indel' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_indel_filtered }
    MERGE_FILTERED_VCFS.out.vcf.filter{ it[0]=='invariant' }.map{ _, vcf, tbi -> [vcf,tbi] }.first().set { ch_invariant_filtered }

    // QC plots for sites and genotypes
    PLOT_VCF_FILTERS (
        FILTER_VCF.out.tables.collect()
    )

    // QC plots for sample missing data
    PLOT_SAMPLE_FILTERS (
        CALC_DATASET_FILTERS.out.missing_summary.first(),
        params.sample_max_missing
    )

    FILTER_VCF.out.samples_to_keep
        .splitText( by: 1 )
        .unique()
        .set { ch_sample_names_filt }

    // Calculate VCF statistics
   VCF_STATS (
        ch_combined_filtered,
        ch_genome_indexed,
        ch_sample_names_filt
    )
        
    // Subset the merged vcf channels to each variant type for emission
    emit:
    filtered_combined = ch_combined_filtered
    filtered_snps = ch_snp_filtered
    filtered_indels = ch_indel_filtered
    reports = VCF_STATS.out.vcfstats

}