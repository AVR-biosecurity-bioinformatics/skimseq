/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype' 
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { CREATE_BED_INTERVALS                                   } from '../modules/create_bed_intervals'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import' 
include { GENOTYPE_POSTERIORS                                    } from '../modules/genotype_posteriors' 


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
    
    // If intervals are being subdivided on the mask, use the mask bed
    if ( params.interval_subdivide_at_masks ){
          ch_breakpoints_bed = ch_mask_bed_gatk
        } else {
          ch_breakpoints_bed = ch_dummy_file
    }
        
    // create groups of genomic intervals for parallel genotyping
    CREATE_BED_INTERVALS (
        ch_genome_indexed,
        ch_include_bed,
        ch_breakpoints_bed,
        params.interval_n,
        params.interval_size,
        params.interval_subdivide_balanced
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_BED_INTERVALS.out.interval_bed
        .flatten()
        // get hash from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_hash = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_bed ] }
        .set { ch_interval_bed }
        

    // combine sample-level bams with each interval_bed file and interval hash
    ch_sample_bam
        .combine ( ch_interval_bed )
        .set { ch_sample_intervals }

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed,
        params.interval_padding,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        params.hc_min_pruning,
        params.hc_min_dangling_length,
        params.hc_max_reads_startpos,
        params.ploidy
    )

    // group GVCFs by interval 
    CALL_VARIANTS.out.gvcf_intervals
        .map { sample, gvcf, tbi, interval_hash, interval_bed -> [ interval_hash, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed, by: 0 )
        .map { interval_hash, gvcf, tbi, interval_bed -> [ interval_hash, interval_bed, gvcf, tbi ] }
        .set { ch_gvcf_interval }

    // Import GVCFs into a genomicsDB per Interval
    GENOMICSDB_IMPORT (
        ch_gvcf_interval,
        ch_genome_indexed
    )

    // joint-call genotypes
    JOINT_GENOTYPE (
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        params.output_invariant
    )

    // calculate genotype posteriors
    GENOTYPE_POSTERIORS (
        JOINT_GENOTYPE.out.vcf
    )

    // collect .vcfs into a single element
    GENOTYPE_POSTERIORS.out.vcf
        .collect()
        .set { ch_vcfs }

    // merge interval .vcfs into a single file
    MERGE_VCFS (
        ch_vcfs
    )

    emit: 
    vcf = MERGE_VCFS.out.vcf
    //posteriors = ch_posteriors


}