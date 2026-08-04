/*
    Genotype samples using GATK
*/

//// import modules
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 
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
        .map { sample, bed, tbi -> tuple(bed, tbi) }   // keep bed+tbi pairs
        .toList()
        .filter { lst -> lst && !lst.isEmpty() }
        .map { pairs ->
            def beds = pairs.collect { it[0] }
            def tbis = pairs.collect { it[1] }
            tuple("joint", beds, tbis)
        }
        .set { ch_counts }

    // Create mpileup intervals
    CREATE_INTERVAL_CHUNKS_MP (
        ch_counts,
        ch_genome_indexed,
        ch_include_bed.first(),
        params.mp_bases_per_chunk,
        params.min_interval_gap,
        params.split_large_intervals,
        "false"
    )

    CREATE_INTERVAL_CHUNKS_MP.out.interval_bed
        .flatMap { sample, beds, tbis  ->
            // normalize to a list for cases where there are only 1 bed output for a sample
            def bedList = (beds instanceof List) ? beds : [beds]
            def tbiList = (tbis instanceof List) ? tbis : [tbis]

            assert bedList.size() == tbiList.size() :
            "Mismatch for ${sample}: beds=${bedList.size()} tbis=${tbiList.size()}"

            // emit one tuple per bed file
            (0..<bedList.size()).collect { i ->
                def bed = bedList[i] as Path
                def tbiPath = tbiList[i]
                def base = bed.getFileName().toString()
                base = base.replaceFirst(/\.gz$/, '')
                base = base.replaceFirst(/\.bed$/, '')
                def interval_hash = base.startsWith('_') ? base.substring(1) : base
                tuple(interval_hash, bed, tbiPath)
            }
        }
        .filter { interval_hash, interval_bed, bed_tbi -> interval_bed && interval_bed.size() > 0 }   // drop empty
        .ifEmpty {
            log.warn(
                "No mpileup intervals remained after coverage & inclusion filtering, " +
                "Variant calling will be skipped. Check your max_depth argument"
            )
            tuple('__NO_INTERVALS__', null, null)
        }
        .filter { interval_hash, interval_bed, bed_tbi ->
            interval_hash != '__NO_INTERVALS__'
        }
        .set { ch_interval_bed_mp }

    // combine sample-level cran with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_cram 
        .combine ( ch_interval_bed_mp )
        .map { sample, cram, crai, interval_chunk, interval_bed, bed_tbi -> [ interval_chunk, cram, crai ] }
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
        ch_popmap.first()
    )
    
    ch_merged_unfiltered_vcf = Channel.empty()
    if ( params.output_unfiltered_vcf ){
        // TODO: Make this output seperate files for each variant type
        MPILEUP.out.vcf
            .map { interval_chunk, interval_bed, bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
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