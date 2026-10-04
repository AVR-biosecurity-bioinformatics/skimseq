/*
    Genotype samples using bcftools
*/

//// import modules
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks/create_interval_chunks'
include { PROJECT_WORKLOAD_TO_INTERVALS                          } from '../modules/project_workload_to_intervals/project_workload_to_intervals'
include { MPILEUP                                                } from '../modules/mpileup/mpileup'

workflow BCFTOOLS_CALLING {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_calling_bed
    ch_popmap
    ch_cohort_size
    ch_included_bases

    main: 

    // Calculate nchunks
    ch_n_chunks = ch_cohort_size
        .combine(ch_included_bases)
        .map { cohort_size, included_bases ->
            long total_bases =
                (cohort_size as long) * (included_bases as long)

            long bases_per_chunk =
                params.mp_target_sample_bases as long

            Math.max(
                1,
                ((total_bases + bases_per_chunk - 1) / bases_per_chunk) as int
            )
        }


    /* 
       Create groups of genomic intervals for parallel calling
       This is done based on CRAI indexes similar to goleft indexsplit
       This is a fast but coarse way of assessign workload
    */

    ch_crai_workload_inputs = ch_sample_cram
        .map { _sample, _cram, crai ->
        crai
        }
        .collect()

    PROJECT_WORKLOAD_TO_INTERVALS(
        ch_crai_workload_inputs,
        ch_calling_bed,
        ch_genome_indexed
    )

    CREATE_INTERVAL_CHUNKS_MP(
        PROJECT_WORKLOAD_TO_INTERVALS.out.bed,
        ch_n_chunks,
        params.split_large_intervals
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
       Joint call variants per chunk
    */

    // Joint calling using mpileup
    MPILEUP (
        ch_cram_interval,
        ch_genome_indexed,
        ch_cohort_size,
        ch_popmap.first()
    )
    
    // Merged unfiltered bcf outputs - just used for publishing
    ch_merged_unfiltered_bcf = channel.empty()
    if ( params.output_unfiltered_bcf ){
        // TODO: Make this output seperate files for each variant type
        MPILEUP.out.bcf
            .map { _interval_chunk, _interval_bed, _bed_tbi, bcf, csi -> tuple('unfiltered', bcf, csi) }
            .groupTuple(by: 0)
            .set { ch_bcf_to_merge }

        CONCAT_UNFILTERED_VCFS (
            ch_bcf_to_merge
        )
    
        CONCAT_UNFILTERED_VCFS.out.vcf
            .set { ch_merged_unfiltered_bcf }
    }

    emit: 
    bcf = MPILEUP.out.bcf
    merged_unfiltered_bcf = ch_merged_unfiltered_bcf

}