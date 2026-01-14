/*
    Filter .vcf files from GATK
*/

//// import modules
include { CALC_CHUNK_MISSING                           } from '../modules/calc_chunk_missing'
include { MERGE_CHUNK_MISSING                          } from '../modules/merge_chunk_missing'
include { FILTER_VCF_GENOTYPES                         } from '../modules/filter_vcf_sites'
include { MERGE_VCFS as MERGE_FILTERED_VCFS            } from '../modules/merge_vcfs'
include { VCF_STATS                                    } from '../modules/vcf_stats'
include { PLOT_VCF_FILTERS                             } from '../modules/plot_vcf_filters'
include { PLOT_SAMPLE_FILTERS                          } from '../modules/plot_sample_filters'

workflow FILTER_GENOTYPES {

    take:
    ch_vcfs
    ch_filtered_sites
    ch_mask_bed_vcf
    ch_sample_names

    main: 


    // For each input VCF, combine with type to make a copy for each variant type, then run FILTER_VCF on each
    def types = ['snp', 'indel']
    if( params.output_invariant ) {
    types << 'invariant'
    }

    // Filter VCF genotypes

    // Calculate per-sample missing data, creating a list of samples to keep


    // Filter genotypes for quality - Set to missing genotype but retain GL/PL for probabilistic analyses
    // Apply final missing data
    FILTER_VCF_GENOTYPES (
        ch_vcfs.combine( channel.of(*types) ),
	    ch_mask_bed_vcf
    )

    // Create a channel of all 3 variant types + all together for merging
    ch_vcfs_nonempty
        .concat(ch_vcfs_nonempty.map { type, vcf, tbi -> tuple('combined', vcf, tbi) })
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

    // Calculate VCF statistics
   VCF_STATS (
        ch_combined_filtered,
        ch_genome_indexed,
        ch_sample_names
    )
        

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