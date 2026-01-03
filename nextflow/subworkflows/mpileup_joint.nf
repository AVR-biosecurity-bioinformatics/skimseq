/*
    Genotype samples using GATK
*/

//// import modules
include { VALIDATE_GVCF                                          } from '../modules/validate_gvcf'
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { COUNT_VCF_DEPTH                                        } from '../modules/count_vcf_depth'
include { CREATE_INTERVAL_CHUNKS_HC                              } from '../modules/create_interval_chunks_hc'
include { PROFILE_HC                                             } from '../modules/profile_hc'
include { STAGE_GVCF                                             } from '../modules/stage_gvcf'

workflow GATK_SINGLE {

    take:
    ch_sample_names
    ch_sample_cram
    ch_rg_to_validate
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed
    ch_dummy_file

    main: 


    /* 
       Create groups of genomic intervals for parallel genotype calling
    */

    // Create haplotypecaller intervals on an all-dataset basis
    CREATE_INTERVAL_CHUNKS_HC (
        ch_sample_cram,
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        params.hc_bases_per_chunk,
        params.split_large_intervals,
        params.hc_rmdup,
        params.hc_minbq,
        params.hc_minmq,
        params.min_interval_gap
    )
   
    // CREATE_INTERVAL_CHUNKS_HC.out.interval_bed emits: tuple(sample, bed)
    // where `bed` is either a List<Path> or a single Path, so has to be normalised to list
    CREATE_INTERVAL_CHUNKS_HC.out.interval_bed
        .flatMap { sample, beds ->
            // normalize to a list for cases where there are only 1 bed output for a sample
            def lst = (beds instanceof List) ? beds : [ beds ]
            // emit one tuple per bed file
            lst.collect { bed ->
            bed  = bed as Path
            def base = bed.baseName
            def interval_chunk = base.startsWith('_') ? base.substring(1) : base
            tuple(sample, interval_chunk, bed)
            }
        }
        .set { ch_interval_bed_hc }

    // Combine intervals with cram files for genotyping
    ch_interval_bed_hc 
	    .combine( ch_sample_cram, by: [0, 0] )
        .set { ch_sample_intervals }

    /* 
       Call variants per sample
    */

    // call variants for single samples across intervals
    MPILEUP (
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_gatk
    )

    // Merge interval GVCFs by sample
    MPILEUP.out.vcf_intervals
        .map { interval_chunk, interval_bed, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
        .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_UNFILTERED_VCFS (
        ch_vcf_to_merge
    )

    emit: 
    gvcf = STAGE_GVCF.out.gvcf
}