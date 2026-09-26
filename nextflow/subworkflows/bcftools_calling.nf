/*
    Genotype samples using GATK
*/

//// import modules
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 
include { COMBINE_MOSDEPTH_EVENTS                                } from '../modules/combine_mosdepth_events/combine_mosdepth_events'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks/create_interval_chunks'
include { MPILEUP                                                } from '../modules/mpileup/mpileup'

workflow BCFTOOLS_CALLING {

    take:
    ch_sample_names
    ch_sample_cram
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_genotype
    ch_read_counts
    ch_popmap

    main: 

    /* 
       Create groups of genomic intervals for parallel genotyping
    */

    ch_read_counts
        .map { _sample, starch -> starch }
        .toList()
        .filter { archives -> !archives.isEmpty() }
        .set { ch_events }

    COMBINE_MOSDEPTH_EVENTS(
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_mask_bed_genotype,
        ch_events,
        "false"
    )

    // cohort_rle.bed column 4 is depth per base.
    // CREATE_INTERVAL_CHUNKS_MP calculates the weight of each emitted span.
    CREATE_INTERVAL_CHUNKS_MP(
        ch_include_bed,
        COMBINE_MOSDEPTH_EVENTS.out.rle,
        params.mp_bases_per_chunk,
        params.min_interval_gap,
        false
    )

    CREATE_INTERVAL_CHUNKS_MP.out.interval_bed
        .flatMap { _name, beds, tbis ->
            def bedList = (beds instanceof List) ? beds : [beds]
            def tbiList = (tbis instanceof List) ? tbis : [tbis]

            // Match indexes by filename, not by position in two glob lists.
            def tbiByName = tbiList.collectEntries { tbi ->
                [(tbi.name): tbi]
            }

            bedList.collect { bed ->
                def tbi = tbiByName["${bed.name}.tbi"]

                assert tbi != null :
                    "Missing tabix index for ${bed.name}"

                def interval_id = bed.name.replaceFirst(/\.bed\.gz$/, '')

                tuple(interval_id, bed, tbi)
            }
        }
        .ifEmpty {
            log.warn(
                "No mpileup intervals remained after coverage filtering; " +
                "variant calling will be skipped."
            )
            tuple('__NO_INTERVALS__', null, null)
        }
        .filter { interval_id, _bed, _tbi ->
            interval_id != '__NO_INTERVALS__'
        }
        .set { ch_interval_bed_mp }

    // combine sample-level cram with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_cram 
        .combine ( ch_interval_bed_mp )
        .map { _sample, cram, crai, interval_chunk, _interval_bed, _bed_tbi -> [ interval_chunk, cram, crai ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_mp, by: 0 )
        .map { interval_chunk, cram, crai, interval_bed, bed_tbi -> [ interval_chunk, interval_bed, bed_tbi, cram, crai ] }
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
        ch_cohort_size,
        ch_popmap.first(),
        ch_mask_bed_genotype
    )
    
    ch_merged_unfiltered_vcf = channel.empty()
    if ( params.output_unfiltered_vcf ){
        // TODO: Make this output seperate files for each variant type
        MPILEUP.out.vcf
            .map { _interval_chunk, _interval_bed, _bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        CONCAT_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )
    
        CONCAT_UNFILTERED_VCFS.out.vcf
            .set { ch_merged_unfiltered_vcf }
    }

    emit: 
    vcf = MPILEUP.out.vcf
    merged_unfiltered_vcf = ch_merged_unfiltered_vcf

}