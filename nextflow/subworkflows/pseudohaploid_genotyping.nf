/*
    Pseudohaploid genotyping using list of provided sites
*/

//// import modules
include { MERGE_VCFS as MERGE_UNFILTERED_VCFS                    } from '../modules/merge_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks'
include { MPILEUP                                                } from '../modules/mpileup'

workflow BCFTOOLS_GENOTYPING {

    take:
    ch_sites_to_genotype
    ch_sample_cram

    main: 

    // Call pseudohaploids per sample
    CALL_PSEUDOHAPLOID (
        ch_sample_cram,
        ch_sites_to_genotype
    )

    // Merge pseudohaploid
    CREATE_INTERVAL_CHUNKS_MP.out.interval_bed
	    .map { sample,interval_bed -> interval_bed }
        .filter { interval_bed -> interval_bed && interval_bed.size() > 0 }   // drop empty
        .collect()
        .flatten()
        // get interval_chunk from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_chunk = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_chunk, interval_bed ] }
        .set { ch_interval_bed_mp }

    // combine sample-level cran with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_cram 
        .combine ( ch_interval_bed_mp )
        .map { sample, cram, crai, interval_chunk, interval_bed -> [ interval_chunk, cram, crai ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_mp, by: 0 )
        .map { interval_chunk, cram, crai, interval_bed -> [ interval_chunk, interval_bed, cram, crai ] }
        .set { ch_cram_interval }

    /* 
       Call variants per sample
    */

    // Calculate cohort size for memory scaling
    ch_cohort_size = ch_sample_names.unique().count()

    // call variants for single samples across intervals
    MPILEUP (
        ch_cram_interval,
        ch_genome_indexed,
        ch_cohort_size
    )

    MPILEUP.out.vcf
        .map { interval_chunk, interval_bed, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
        .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_UNFILTERED_VCFS (
        ch_vcf_to_merge
    )

    // How to set up count_vcf_records in a way thats flexible to mpileup as well?
    
    emit: 
    vcf = MPILEUP.out.vcf
    merged_vcf = MERGE_UNFILTERED_VCFS.out.vcf
    //missing_frac = COUNT_VCF_RECORDS.out.missing_frac
    //variant_dp = COUNT_VCF_RECORDS.out.variant_dp
}