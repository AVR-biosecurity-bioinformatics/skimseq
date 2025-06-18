/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                                 } from '../modules/call_variants'
include { COMBINE_GVCFS                                              } from '../modules/combine_gvcfs' 
include { CONVERT_INTERVALS                                             } from '../modules/convert_intervals' 
include { CREATE_BEAGLE                                              } from '../modules/create_beagle' 
include { CREATE_INTERVALS                                              } from '../modules/create_intervals' 
include { GENOTYPE_POSTERIORS                                              } from '../modules/genotype_posteriors' 
include { JOINT_GENOTYPE                                              } from '../modules/joint_genotype' 
include { MERGE_VCFS                                              } from '../modules/merge_vcfs' 



workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed

    main: 

    // check that params.interval_size is not less than 100k bases
    // Disabled - this could be covered in documentation
    //if ( params.interval_size.toFloat() < 1E5 ){
    //    println ("\n*** ERROR: 'params.interval_size' must be >= 1E5 (100,000 bases) ***\n")
    //    error()
    //}

    /* 
        Genotype samples individually and jointly
    */

    // create genome intervals for genotyping
    CREATE_INTERVALS (
        ch_genome_indexed,
        params.interval_size
    )

    // split intervals file into chunks of 50 lines for conversion
    CREATE_INTERVALS.out.intervals
        .splitText ( by: 50, file: true )
        .set { ch_intervals }

    // turn intervals into GATK format via .bed
    CONVERT_INTERVALS (
        ch_intervals,
        ch_genome_indexed
    )

    // create intervals channel, with one interval_list file per element
    CONVERT_INTERVALS.out.interval_list
        .flatten()
        // get hash from interval_list name as element to identify intervals
        .map { interval_list ->
            def interval_hash = interval_list.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_list ] }
        .set { ch_interval_list }

    // combine sample-level bams with each interval_list file and interval hash
    ch_sample_bam
        .combine ( ch_interval_list )
        .set { ch_sample_intervals }

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed
    )

    // group GVCFs by interval 
    CALL_VARIANTS.out.gvcf_intervals
        .map { sample, gvcf, tbi, interval_hash, interval_list -> [ interval_hash, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_list, by: 0 )
        .map { interval_hash, gvcf, tbi, interval_list -> [ interval_hash, interval_list, gvcf, tbi ] }
        .set { ch_gvcf_interval }

    // combine GVCFs into one file per interval
    COMBINE_GVCFS (
        ch_gvcf_interval,
        ch_genome_indexed
    )

    // calculate genotype posteriors over each genomic interval
    GENOTYPE_POSTERIORS (
        COMBINE_GVCFS.out.gvcf_intervals
    )

    // call genotypes at variant sites
    JOINT_GENOTYPE (
        GENOTYPE_POSTERIORS.out.gvcf_intervals,
        ch_genome_indexed
    )

    // collect .vcfs into a single element
    JOINT_GENOTYPE.out.vcf
        .collect()
        .set { ch_vcfs }

    // merge interval .vcfs into a single file
    MERGE_VCFS (
        ch_vcfs
    )

    // get just posterior .g.vcfs from channel
    GENOTYPE_POSTERIORS.out.gvcf_intervals
        .map { interval_hash, interval_list, gvcf, gvcf_tbi -> [ gvcf, gvcf_tbi ] }
        .set { ch_posteriors }

    // create beagle file from .g.vcf files with posteriors
    /// NOTE: Could move this to an ANGSD-specific workflow
    CREATE_BEAGLE (
        ch_posteriors,
        ch_genome_indexed
    )

    emit: 
    vcf = MERGE_VCFS.out.vcf
    posteriors = ch_posteriors


}