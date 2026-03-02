/*
    Genotype samples using GATK
*/

//// import modules
include { MERGE_VCFS as MERGE_UNFILTERED_VCFS                    } from '../modules/merge_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks'
include { MPILEUP                                                } from '../modules/mpileup'

workflow MPILEUP_CALLING {

    take:
    ch_sample_names
    ch_sample_cram
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_genotype
    ch_read_counts

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

    // Get Mpileup parameters
    def MPILEUP_DEFAULTS = [
        RMDUP:'NA', EXCLUDE_PAD:'NA', OUTPUT_INVARIANT:'NA', PLOIDY:'NA',
        MINBQ:'NA', MINMQ:'NA', MIN_ALIGNED_LENGTH:'NA',
        MIN_FRAGMENT_LENGTH:'NA', MAX_FRAGMENT_LENGTH:'NA',
        MUTATION_RATE:'NA', MAXDEPTH:'NA'
    ]

    def canon = { Map m, Map defaults = [:] ->
        def merged = defaults + m
        merged.collect { k,v -> "${k}=${v == null ? 'NA' : v}" }.sort().join(';')
    }

    def mpileupMap = [
        RMDUP: params.rmdup,
        EXCLUDE_PAD: params.exclude_padding,
        OUTPUT_INVARIANT: params.output_invariant,
        PLOIDY: params.ploidy,
        MINBQ: params.minbq,
        MINMQ: params.minmq,
        MIN_ALIGNED_LENGTH: params.min_aligned_length,
        MIN_FRAGMENT_LENGTH: params.min_fragment_length,
        MAX_FRAGMENT_LENGTH: params.max_fragment_length,
        MUTATION_RATE: params.mutation_rate,
        MAXDEPTH: params.max_depth
    ]

    ch_mpileup_kv = Channel.value( canon(mpileupMap, MPILEUP_DEFAULTS) )

    // call variants for single samples across intervals
    MPILEUP (
        ch_cram_interval,
        ch_genome_indexed,
        ch_cohort_size,
        ch_mpileup_kv
    )
    
    if ( params.output_unfiltered_vcf ){
        // TODO: Make this output seperate files for each variant type
        MPILEUP.out.vcf
            .map { interval_chunk, interval_bed, bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        MERGE_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )
    }

    emit: 
    vcf = MPILEUP.out.vcf

}