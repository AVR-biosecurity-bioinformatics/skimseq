/*
    Genotype samples using GATK
*/

//// import modules
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype' 
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { MERGE_VCFS as MERGE_UNFILTERED_VCFS                    } from '../modules/merge_vcfs' 
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_SHORT           } from '../modules/count_vcf_records'
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_LONG            } from '../modules/count_vcf_records'
include { COUNT_VCF_DEPTH                                        } from '../modules/count_vcf_depth'
include { CREATE_INTERVAL_CHUNKS_JC                              } from '../modules/create_interval_chunks_jc'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import' 
include { PROFILE_JC                                             } from '../modules/profile_jc' 

workflow GATK_JOINT {

    take:
    ch_sample_gvcf
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed
    ch_dummy_file
    ch_sample_names

    main: 

    // Count missing data in each gvcf - this will be used later for missing data and percentile depth filtering
    COUNT_VCF_DEPTH (
        ch_sample_gvcf,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )

    /* 
       Create groups of genomic intervals for parallel joint calling
    */

    // Count number of reads contained within each interval, for long contigs (chromosomes)
    // Note: For the long contigs use the mask to split them into smaller genotypable chunks
    COUNT_VCF_RECORDS_LONG (
        ch_sample_gvcf,
        ch_long_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .set { counts_long }

    // Count number of vcf records contained within each interval, for short contigs (scaffolds)
    // NOTE: For the short contigs use the full contigs with NO MASK, by providing ch_dummy_file
    // This ensures compatibility with --merge-contigs-into-num-partitions in genomicsdbimport
    // Masked regions will still not be included if the mask was used for call_variants
    COUNT_VCF_RECORDS_SHORT (
        ch_sample_gvcf,
        ch_short_bed.first(),
        ch_dummy_file,
        ch_genome_indexed
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .set { counts_short }

    // Merge both long and short bed list channels
    counts_long
        .mix(counts_short)
        .filter { L -> L && L.size() > 0 } // drop empty [] lists
        .set { ch_long_short_beds }

    // Create joint calling intervals, long and short processed separately
    // Takes the sum of counts * samples - i.e. number of genotypes
    CREATE_INTERVAL_CHUNKS_JC (
        ch_long_short_beds,
        params.jc_genotypes_per_chunk
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
    ch_sample_gvcf 
        .combine ( ch_interval_bed_jc )
        .map { sample, gvcf, tbi, interval_hash, interval_bed -> [ interval_hash, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_jc, by: 0 )
        .map { interval_hash, gvcf, tbi, interval_bed -> [ interval_hash, interval_bed, gvcf, tbi ] }
        .set { ch_gvcf_interval }


    // Calculate cohort size from sample names - This is ued for memory scaling of next 2 steps
    ch_cohort_size = ch_sample_names.unique().count()

    // Import GVCFs into a genomicsDB per Interval
    GENOMICSDB_IMPORT (
        ch_gvcf_interval,
        ch_genome_indexed,
        ch_cohort_size
    )

    // joint-call genotypes across all samples per Interval
    JOINT_GENOTYPE (
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        params.output_invariant,
        ch_cohort_size
    )

    if( params.profile_gatk ) {
        // TODO: if extra covariates are calculated for these intervals, run separately and then merge output TSV in nextflow
        PROFILE_JC (
            JOINT_GENOTYPE.out.log.map { interval_hash, log -> log}.collect()
        )
    }

    if( params.output_unfiltered_vcf ) {
        JOINT_GENOTYPE.out.vcf
            .map { interval_hash, interval_bed, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        MERGE_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )
        
    }

    emit: 
    vcf = JOINT_GENOTYPE.out.vcf
    missing_frac = COUNT_VCF_DEPTH.out.missing_frac
    variant_dp = COUNT_VCF_DEPTH.out.variant_dp
}