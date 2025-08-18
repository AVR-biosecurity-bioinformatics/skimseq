/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype' 
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { CREATE_HC_INTERVALS                                    } from '../modules/create_bed_intervals'
include { CREATE_JC_INTERVALS                                    } from '../modules/create_jc_intervals'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import' 

workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk

    main: 

    /* 
        Genotype samples individually and jointly
    */
            
    // create groups of genomic intervals for parallel haplotypecaller
    // TODO: Replace this with a process that takes into account the number of aligned bases
    CREATE_BED_INTERVALS_HC (
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_gatk,
        params.hc_interval_n,
        params.hc_interval_size,
        params.interval_subdivide_balanced
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_BED_INTERVALS_HC.out.interval_bed
        .flatten()
        // get hash from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_hash = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_bed ] }
        .set { ch_interval_bed_hc }
        

    // combine sample-level bams with each interval_bed file and interval hash
    ch_sample_bam
        .combine ( ch_interval_bed_hc )
        .set { ch_sample_intervals }

    // collect haplotypecaller parameters into a single list
    Channel.of(
        params.hc_min_pruning,
        params.hc_min_dangling_length,
        params.hc_max_reads_startpos,
        params.hc_rmdup,
        params.hc_minmq,
        params.ploidy
    )
    .collect( sort: false )
    .set { ch_hc_params }

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed,
        params.interval_padding,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        ch_hc_params
    )

    // Merge interval GVCFs by sample
    // TODO: This would be faster using GatherGVCFs, but need to make sure the intervals are ordered
    CALL_VARIANTS.out.gvcf_intervals
        .groupTuple ( by: 0 )
        .set { ch_gvcf_to_merge }

    MERGE_GVCFS (
        ch_gvcf_to_merge.map { sample, gvcf, tbi, interval_hash, interval_bed -> [ gvcf, tbi ] },
        ch_gvcf_to_merge.flatMap { sample, gvcf, tbi, interval_hash, interval_bed -> [ sample ] }
    )

    // Calculate number of intervals to make for parallel joint calling
    /* number of chunks = round( params.jc_interval_sample_scale * N ), clamped to ≥ 2 */
    MERGE_GVCFS.out.vcf
    .count()
    .map { n ->
        double f = (params.jc_interval_sample_scale as double)
        assert f > 0 && f <= 1 : "params.chunk_frac must be in (0,1], got: ${f}"
        int nchunks = Math.max(2, (int)Math.round(n * f))
        nchunks
    }
    .set { ch_jc_nchunks }   // => emits a single integer

    // create groups of genomic intervals for parallel joint calling
    CREATE_JC_INTERVALS (
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_gatk,
        ch_jc_nchunks
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_JC_INTERVALS.out.interval_bed
        .flatten()
        // get hash from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_hash = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_bed ] }
        .set { ch_interval_bed_jc }
        

    // combine sample-level gvcf with each interval_bed file and interval hash
    // Then group by interval for joint genotyping
    MERGE_GVCFS.out.vcf 
        .combine ( ch_interval_bed_jc )
        .map { gvcf, tbi, interval_hash, interval_bed -> [ interval_hash, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_jc, by: 0 )
        .map { interval_hash, gvcf, tbi, interval_bed -> [ interval_hash, interval_bed, gvcf, tbi ] }
        .set { ch_gvcf_interval }

    // Import GVCFs into a genomicsDB per Interval
    GENOMICSDB_IMPORT (
        ch_gvcf_interval,
        ch_genome_indexed
    )

    // joint-call genotypes across all samples per Interval
    JOINT_GENOTYPE (
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        params.output_invariant
    )

    // collect joint called .vcfs into a single element
    // Use collect without flattening and transpose to get [[vcf] [tbi]]
    JOINT_GENOTYPE.out.vcf
        .map { interval_hash, interval_bed, vcf, tbi -> [ vcf, tbi ] }
        .collect(flat: false)
 	    .map { it.transpose() }
        .set { ch_vcfs }

    // merge interval .vcfs into a single file
    MERGE_VCFS (
        ch_vcfs,
        "joint"
    )

    emit: 
    vcf = MERGE_VCFS.out.vcf
    //posteriors = ch_posteriors


}