/*
    Pseudohaploid genotyping using list of provided sites
*/

//// import modules
include { MERGE_VCFS as MERGE_SITELIST                           } from '../modules/merge_vcfs' 
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
        .map { sample, cram, crai, variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, sites_vcf, sites_tbi -> [ variant_type, interval_hash, cram, crai ] }
        .groupTuple ( by: [0,1] )
        // join to get back interval_file
        .join ( ch_sites_to_genotype, by: [0,1] )
        // variant type and interval hash columns are combined into a single string for compatibility with mpileup
        .map { variant_type, interval_hash, cram, crai, interval_bed, bed_tbi, vcf, tbi, sites_vcf, sites_tbi -> 
            def id = "${variant_type}_${interval_hash}"
            tuple(id, sites_vcf, sites_tbi, cram, crai) 
        }
	    .set { ch_cram_to_genotype }

    // Calculate cohort size for memory scaling
    ch_cohort_size = ch_sample_names.unique().count()

    // Call just target sites using mpileup
    MPILEUP_PSEUDOHAP (
        ch_cram_to_genotype,
        ch_genome_indexed,
        ch_cohort_size
    )
     
    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
            MPILEUP_PSEUDOHAP.out.vcf.map { id, sites_vcf, sites_tbi, vcf, tbi -> tuple(id, vcf, tbi) },
            ch_genome_indexed
    )

        

    // Split the variant_type and interval_hash back out to separate columns
    CREATE_PSEUDOHAP.out.vcf
    .map { id, vcf, tbi ->
        def (variant_type, interval_hash) = id.split('_', 2)
        tuple(variant_type, interval_hash, vcf, tbi)
    }
    .join ( ch_sites_to_genotype.map {
         variant_type, interval_hash, interval_bed, bed_tbi, vcf, tbi, sites_vcf, sites_tbi
          -> tuple(variant_type, interval_hash, sites_vcf, sites_tbi)
        }, by: [0,1] ) 
    .map { variant_type, interval_hash, vcf, tbi, sites_vcf, sites_tbi -> tuple(variant_type, interval_hash, sites_vcf, sites_tbi, vcf, tbi) }
    .set{ ch_pseudohaploid_vcf }

    emit: 
    vcf = ch_pseudohaploid_vcf
}