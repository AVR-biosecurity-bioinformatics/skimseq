/*
    Pseudohaploid genotyping using list of provided sites
*/

//// import modules
include { MERGE_VCFS as MERGE_PSEUDOHAP                          } from '../modules/merge_vcfs' 
include { MPILEUP as MPILEUP_PSEUDOHAP                           } from '../modules/mpileup'
include { CREATE_PSEUDOHAP                                       } from '../modules/create_pseudohap' 

workflow PSEUDOHAPLOID_GENOTYPING {

    take:
    ch_sites_to_genotype
    ch_sample_cram
    ch_genome_indexed
    ch_sample_names

    main: 

    // combine sample-level cram with each interval_bed file and interval chunk
    ch_sample_cram 
        .combine ( ch_sites_to_genotype )
        .map { sample, cram, crai, interval_hash, interval_bed, bed_tbi, vcf, tbi, sites_vcf, sites_tbi -> [ interval_hash, cram, crai ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_sites_to_genotype, by: 0 )
        .map { interval_hash, cram, crai, interval_bed, bed_tbi, vcf, tbi, sites_vcf, sites_tbi -> [ interval_hash, sites_vcf, sites_tbi, cram, crai ] }
	    .set { ch_cram_to_genotype }

    // Calculate cohort size for memory scaling
    ch_cohort_size = ch_sample_names.unique().count()

    // Call just target sites using mpileup
    MPILEUP_PSEUDOHAP (
        ch_cram_to_genotype,
        ch_genome_indexed,
        ch_cohort_size
    )
    
    MPILEUP_PSEUDOHAP.out.vcf
        .map { interval_hash, sites_vcf, sites_tbi, vcf, vcf_tbi -> [ interval_hash, vcf, vcf_tbi ] }
        .set { ch_vcf_for_pseudo }


    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
            ch_vcf_for_pseudo,
            ch_genome_indexed
    )

    CREATE_PSEUDOHAP.out.vcf
        .map { interval_hash, vcf, tbi -> tuple('pseudohaploid', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_PSEUDOHAP (
        ch_vcf_to_merge
    )

    emit: 
    vcf = MERGE_PSEUDOHAP.out.vcf
}