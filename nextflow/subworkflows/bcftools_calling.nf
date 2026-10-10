nextflow.enable.types = true

/*
    Genotype samples using bcftools
*/


//// import modules
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks/create_interval_chunks'
include { PROJECT_WORKLOAD_TO_INTERVALS                          } from '../modules/project_workload_to_intervals/project_workload_to_intervals'
include { MPILEUP                                                } from '../modules/mpileup/mpileup'

// Import record types
include {
    SampleCram
    PreparedGenome
    VariantChunk
    VariantFile
} from '../types'

workflow BCFTOOLS_CALLING {

    take:
    sample_cram: Channel<SampleCram>
    genome: Value<PreparedGenome>
    popmap: Value<Path>
    cohort_size: Value<Integer>

    main: 
    included_bases = genome
        .map { prepared ->
            prepared.reference_bases.text.trim().toLong()
        }

    // Calculate number of chunks for parallel calling
    n_chunks = cohort_size
        .combine(included_bases)
        .map { cohort_size_value, included_bases_value ->

            def total_bases =
                cohort_size_value * included_bases_value

            def bases_per_chunk =
                params.mp_target_sample_bases

            Math.max(
                1,
                ((total_bases + bases_per_chunk - 1) / bases_per_chunk) as int
            )
        }

    /*
    * Estimate cohort workload from the CRAM indexes.
    */

    crai_workload_inputs = sample_cram
        .map { sample ->
            sample.crai
        }
        .collect()

    projected_workload = PROJECT_WORKLOAD_TO_INTERVALS(
        crai_workload_inputs,
        genome
    )

    interval_chunks = CREATE_INTERVAL_CHUNKS_MP(
        projected_workload,
        n_chunks,
        params.split_large_intervals
    )

    /*
    * Expand the chunk result into one record per interval.
    *
    * CREATE_INTERVAL_CHUNKS returns:
    *
    * record(
    *     name: String,
    *     beds: Set<Path>,
    *     tbis: Set<Path>
    * )
    */

    interval_beds = interval_chunks
        .flatMap { chunks ->
            chunks.beds.collect { bed ->

                def tbi = chunks.tbis.find { candidate ->
                    candidate.name == "${bed.name}.tbi"
                }

                if (tbi == null) {
                    throw new IllegalStateException(
                        "Missing tabix index for ${bed.name}"
                    )
                }

                record(
                    interval_id: bed.name.replaceFirst(/\.bed\.gz$/, ''),
                    interval_bed: bed,
                    bed_tbi: tbi
                )
            }
        }

    /*
    * Form the Cartesian product:
    *
    *     every sample CRAM × every interval
    *
    * combine() returns a tuple containing:
    *
    *     SampleCram, IntervalRecord
    */

    sample_interval_pairs = sample_cram
        .combine(interval_beds)
        .map { sample, interval ->
            tuple(
                interval.interval_id,
                sample
            )
        }

    samples_by_interval = sample_interval_pairs
        .groupBy()

    interval_samples = samples_by_interval
        .map { interval_id, samples ->
            record(
                interval_id: interval_id,
                samples: samples
            )
        }

    interval_records = interval_beds
        .map { interval ->
            record(
                interval_id: interval.interval_id,
                interval_bed: interval.interval_bed,
                interval_index: interval.bed_tbi
            )
        }

    cram_intervals = interval_samples
        .join(
            interval_records,
            by: 'interval_id'
        )
        .map { result ->
            def samples = result.samples
                .toSorted { sample ->
                    sample.sample
                }

            record(
                interval_id: result.interval_id,
                interval_bed: result.interval_bed,
                interval_index: result.interval_index,
                crams: samples.collect { sample ->
                    sample.cram
                },
                crais: samples.collect { sample ->
                    sample.crai
                }
            )
        }
    /*
    * Joint variant calling, once per interval.
    */

    variant_chunks = MPILEUP(
        cram_intervals,
        genome,
        cohort_size,
        popmap
    )

    /*
    * Optionally concatenate unfiltered variants for publication.
    */

    merged_unfiltered_variants = channel.empty()

    if (params.output_unfiltered_bcf) {

        variants_to_merge = variant_chunks
            .map { result ->
                tuple(
                    'unfiltered',
                    record(
                        variant: result.variants,
                        index: result.variants_index
                    )
                )
            }
            .groupBy()
            .map { outname, files ->
                record(
                    outname: outname,
                    variants: files.collect { result ->
                        result.variant
                    },
                    indices: files.collect { result ->
                        result.index
                    }
                )
            }

        merged_unfiltered_variants = CONCAT_UNFILTERED_VCFS(
            variants_to_merge
        )
    }

    /*
    * Outputs.
    */

    emit:
    variants: Channel<VariantChunk> = variant_chunks
    merged_unfiltered_variants: Channel<VariantFile> = merged_unfiltered_variants

}