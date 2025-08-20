/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype' 
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { COUNT_READS_BED                                        } from '../modules/count_reads_bed'
include { COUNT_VCF_BED as COUNT_VCF_BED_SHORT                   } from '../modules/count_vcf_bed'
include { COUNT_VCF_BED as COUNT_VCF_BED_LONG                    } from '../modules/count_vcf_bed'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_HC    } from '../modules/create_interval_chunks'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_JC    } from '../modules/create_interval_chunks'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import' 

workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed

    main: 

    /* 
       Create groups of genomic intervals for parallel haplotypecaller
    */

    // Count number of reads in each interval
    COUNT_READS_BED (
        ch_sample_bam,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )

    // Create haplotypecaller intervals
    // As next step is parallelised by sample*interval, use the 'mean' mode
    CREATE_INTERVAL_CHUNKS_HC (
        params.hc_bases_per_chunk,
        COUNT_READS_BED.out.counts.map { sample, counts -> [ counts ] }.collect(),
        "mean"
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_INTERVAL_CHUNKS_HC.out.interval_bed
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

    // Merge interval GVCFs by sample
    // TODO: This would be faster using GatherGVCFs, but need to make sure the intervals are ordered
    CALL_VARIANTS.out.gvcf_intervals
        .groupTuple ( by: 0 )
        .set { ch_gvcf_to_merge }

    MERGE_GVCFS (
        ch_gvcf_to_merge.map { sample, gvcf, tbi, interval_hash, interval_bed -> [ gvcf, tbi ] },
        ch_gvcf_to_merge.flatMap { sample, gvcf, tbi, interval_hash, interval_bed -> [ sample ] }
    )

    /* 
       Create groups of genomic intervals for parallel joint calling
    */

    // Count number of reads contained within each interval, for long contigs (chromosomes)
    COUNT_VCF_BED_LONG (
        MERGE_GVCFS.out.vcf,
        ch_long_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )

    // Count number of reads contained within each interval, for short contigs (scaffolds)
    COUNT_VCF_BED_SHORT (
        MERGE_GVCFS.out.vcf,
        ch_short_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )

    COUNT_VCF_BED.out.counts
        .map { sample, counts -> [ counts ] }
        .collect() 
        .concat( 
            COUNT_VCF_BED.out.counts
                .map { sample, counts -> [ counts ] }
                .collect()
        ).set { ch_long_short_beds }


    // Create joint calling intervals, run one job for short 
    // As next step is parallelised by interval, use the 'sum' mode
    CREATE_INTERVAL_CHUNKS_JC (
        params.jc_genotypes_per_chunk,
        ch_long_short_beds,
        "sum"
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_INTERVAL_CHUNKS_JC.out.interval_bed
        .collect()
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
        .map { sample, gvcf, tbi, interval_hash, interval_bed -> [ interval_hash, gvcf, tbi ] }
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
    vcf = MERGE_VCFS.out.vcf.map { sample, vcf, tbi -> [ vcf, tbi ] }
}