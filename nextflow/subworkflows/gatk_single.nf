/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { COUNT_BAM_READS                                        } from '../modules/count_bam_reads'
include { COUNT_VCF_DEPTH                                        } from '../modules/count_vcf_depth'
include { CREATE_INTERVAL_CHUNKS_HC                              } from '../modules/create_interval_chunks_hc'
include { PROFILE_HC                                             } from '../modules/profile_hc'

workflow GATK_SINGLE {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed
    ch_dummy_file

    main: 

    /* 
       Create groups of genomic intervals for parallel haplotypecaller
    */

    // Count number of reads in each interval
    COUNT_BAM_READS (
        ch_sample_cram,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed,
        "bases",
        params.hc_rmdup,
        params.hc_minmq

    )

    // Create haplotypecaller intervals on per sample basis
    CREATE_INTERVAL_CHUNKS_HC (
        COUNT_BAM_READS.out.counts,
        params.hc_bases_per_chunk
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

    // collect haplotypecaller parameters into a single list
    Channel.of(
        params.hc_interval_padding,
        params.hc_min_pruning,
        params.hc_min_dangling_length,
        params.hc_max_reads_startpos,
        params.hc_rmdup,
        params.hc_minbq,
        params.hc_minmq,
        params.ploidy,
    )
    .collect( sort: false )
    .set { ch_hc_params }

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        ch_hc_params
    )

    if( params.profile_gatk ) {

        // Join back onto cram and gvcf based on first 3 columns
        CALL_VARIANTS.out.log
            .join( ch_sample_intervals.map { sample, interval_chunk, interval_bed,cram, crai -> tuple(sample, interval_chunk, cram, crai) }, by:[0,1] )
            .join( CALL_VARIANTS.out.gvcf_intervals, by:[0,1] )
            .map { sample, interval_chunk, logfile, assembly_regions, cram, crai, gvcf, tbi -> tuple(sample, interval_chunk, cram, crai, gvcf, tbi, logfile, assembly_regions ) }
            .set { ch_for_profile }

        // Profile HC runtimes per Sample x Interval
        PROFILE_HC (
            ch_for_profile,
            ch_genome_indexed
        )

        // Merge and output HC profiles
        PROFILE_HC.out.summary
            .collectFile(
                name: 'hc_profiles.tsv',
                storeDir: "${launchDir}/output/gatk_profiles",
                skip: 1,
                keepHeader: true,
                newLine: false,
                sort: true
            )
    }

    // Merge interval GVCFs by sample
    // TODO: This would be faster using GatherGVCFs, but need to make sure the intervals are ordered
    CALL_VARIANTS.out.gvcf_intervals
        .groupTuple ( by: 0 )
        .set { ch_gvcf_to_merge }

    MERGE_GVCFS (
        ch_gvcf_to_merge.map { sample, interval_chunk, gvcf, tbi -> [ sample, gvcf, tbi ] }
    )

    emit: 
    gvcf = MERGE_GVCFS.out.vcf
}