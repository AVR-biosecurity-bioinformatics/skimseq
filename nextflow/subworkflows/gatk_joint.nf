/*
    Genotype samples using GATK
*/

//// import modules
include { JOINT_GENOTYPE                                                 } from '../modules/joint_genotype' 
include { MERGE_VCFS as MERGE_GVCFS                                      } from '../modules/merge_vcfs' 
include { MERGE_VCFS as MERGE_UNFILTERED_VCFS                            } from '../modules/merge_vcfs' 
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_SHORT                   } from '../modules/count_vcf_records'
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_LONG                    } from '../modules/count_vcf_records'
include { CREATE_INTERVAL_CHUNKS_JC as CREATE_INTERVAL_CHUNKS_JC_LONG    } from '../modules/create_interval_chunks_jc'
include { CREATE_INTERVAL_CHUNKS_JC as CREATE_INTERVAL_CHUNKS_JC_SHORT   } from '../modules/create_interval_chunks_jc'
include { GENOMICSDB_IMPORT                                              } from '../modules/genomicsdb_import' 
include { PROFILE_JC                                                     } from '../modules/profile_jc' 

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

    /* 
       Create groups of genomic intervals for parallel joint calling
    */

    // Count number of VCF records contained within each interval, for long contigs (chromosomes)
    // Note: For the long contigs use the mask and min_interval_gap to ensure there are many smaller intervals which will be assigned to parallel chunks in next step
    COUNT_VCF_RECORDS_LONG (
        ch_sample_gvcf,
        ch_long_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed,
        params.min_interval_gap
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .filter { lst -> lst && !lst.isEmpty() }     // prevent empty-list emission
    .set { counts_long }

    // Create joint calling intervals for long beds
    // Takes the sum of vcf records * samples - i.e. number of genotypes to assign intervals to parallel chunks
    // NOTE: split_large_intervals is used here to allow further splitting of intervals that are over params.jc_genotypes_per_chunk
    CREATE_INTERVAL_CHUNKS_JC_LONG (
        counts_long,
        params.jc_genotypes_per_chunk,
        params.split_large_intervals
    )

    // Count number of vcf records contained within each interval, for short contigs (scaffolds)
    // NOTE: For short use the full contig length without splitting, by providing ch_dummy_file as mask, and min_chr_length as the min interval split size
    // This ensures compatibility with --merge-contigs-into-num-partitions in genomicsdbimport
    COUNT_VCF_RECORDS_SHORT (
        ch_sample_gvcf,
        ch_short_bed.first(),
        ch_dummy_file,
        ch_genome_indexed,
        params.min_chr_length
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .filter { lst -> lst && !lst.isEmpty() }     // prevent empty-list emission
    .set { counts_short }

    // Create joint calling intervals for short chunks
    // Takes the sum of vcf records * samples - i.e. number of genotypes to assign intervals to parallel chunks
    // NOTE: set split_large_intervals to FALSE here as --merge-contigs-into-num-partitions 1 requires full contigs
    CREATE_INTERVAL_CHUNKS_JC_SHORT (
        counts_short,
        params.jc_genotypes_per_chunk,
        "false"
    )

    // create intervals channel, with one interval_bed file per element
    // Mix the long contig chunk channels with the short ones - split long and whole short contigs should never be together
    CREATE_INTERVAL_CHUNKS_JC_LONG.out.interval_bed
        .mix(CREATE_INTERVAL_CHUNKS_JC_SHORT.out.interval_bed)
        .collect()
        .flatten()
        // get interval_chunk from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_chunk = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_chunk, interval_bed ] }
        .set { ch_interval_bed_jc }

    // combine sample-level gvcf with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_gvcf 
        .combine ( ch_interval_bed_jc )
        .map { sample, gvcf, tbi, interval_chunk, interval_bed -> [ interval_chunk, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_jc, by: 0 )
        .map { interval_chunk, gvcf, tbi, interval_bed -> [ interval_chunk, interval_bed, gvcf, tbi ] }
        .set { ch_gvcf_interval }


    // Calculate cohort size from sample names
    // NOTE: This is used for memory scaling of GENOMICSDB_IMPORT and JOINT_GENOTYPE which are primarily driven by sample size
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
        // Profile JC runtimes per interval
        PROFILE_JC (
            JOINT_GENOTYPE.out.log,
            ch_genome_indexed
        )        
        // Merge and output JC profiles
        PROFILE_JC.out.summary
            .collectFile(
                name: 'jc_profiles.tsv',
                storeDir: "${launchDir}/output/gatk_profiles",
                skip: 1,
                keepHeader: true,
                newLine: false,
                sort: true
            )
    }

    JOINT_GENOTYPE.out.vcf
        .map { interval_chunk, interval_bed, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
        .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_UNFILTERED_VCFS (
        ch_vcf_to_merge
    )

    emit: 
    vcf = JOINT_GENOTYPE.out.vcf
    merged_vcf = MERGE_UNFILTERED_VCFS.out.vcf
}