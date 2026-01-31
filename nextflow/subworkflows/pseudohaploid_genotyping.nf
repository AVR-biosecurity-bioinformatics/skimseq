/*
    Pseudohaploid genotyping using list of provided sites
*/

//// import modules
include { MERGE_VCFS as MERGE_PSEUDOHAP                          } from '../modules/merge_vcfs' 
include { MPILEUP as MPILEUP_PSEUDOHAP                           } from '../modules/mpileup'
include { CREATE_PSEUDOHAP                                       } from '../modules/create_pseudohap' 

workflow BCFTOOLS_GENOTYPING {

    take:
    ch_sites_to_genotype
    ch_sample_cram

    main: 

    // combine sample-level cran with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_cram 
        .combine ( ch_sites_to_genotype )
        .map { sample, cram, crai, interval_hash, interval_bed, vcf, tbi, sites_vcf, sites_tbi -> [ interval_hash, sites_vcf, sites_tbi, cram, cram_index ] }
        .set { ch_cram_to_genotype}

    // Call just target sites using mpileup
    MPILEUP_PSEUDOHAP (
        ch_cram_to_genotype,
        ch_genome_indexed,
        ch_cohort_size
    )

    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
            MPILEUP.out.vcf,
            ch_genome_indexed
    )

    CREATE_PSEUDOHAP.out.vcf
        .map { interval_chunk, interval_bed, bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_PSEUDOHAP (
        ch_vcf_to_merge
    )

    emit: 
    vcf = MERGE_PSEUDOHAP.out.vcf
}