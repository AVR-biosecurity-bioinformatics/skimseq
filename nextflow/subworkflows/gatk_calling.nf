nextflow.enable.types = true

/*
    Genotype samples using GATK
*/

//// import modules
include { VALIDATE_GVCF                                          } from '../modules/validate_gvcf/validate_gvcf'
include { PROJECT_WORKLOAD_TO_INTERVALS as PROJECT_LONG          } from '../modules/project_workload_to_intervals/project_workload_to_intervals'
include { PROJECT_WORKLOAD_TO_INTERVALS as PROJECT_SHORT         } from '../modules/project_workload_to_intervals/project_workload_to_intervals'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_HC    } from '../modules/create_interval_chunks/create_interval_chunks'
include { COMBINE_BEDS as COMBINE_WORKLOAD_BEDS                  } from '../modules/combine_beds/combine_beds' 
include { SUBSET_BED_TO_INTERVALS                                } from '../modules/subset_bed_to_intervals/subset_bed_to_intervals' 
include { CREATE_JC_BED_FROM_HC                                  } from '../modules/create_jc_bed_from_hc/create_jc_bed_from_hc' 
include { HAPLOTYPECALLER                                        } from '../modules/haplotypecaller/haplotypecaller'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import/genomicsdb_import' 
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype/joint_genotype' 
//include { STAGE_GVCF                                             } from '../modules/stage_gvcf/stage_gvcf'
include { CONCAT_VCFS as CONCAT_GVCFS_JC                         } from '../modules/concat_vcfs/concat_vcfs'
include { CONCAT_VCFS as CONCAT_GVCFS_SAMPLE                     } from '../modules/concat_vcfs/concat_vcfs'
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 

include {
    SampleReads
    SampleCram
    PreparedGenome
    VariantChunk
    VariantFile
    SampleGvcf 
} from '../types'


workflow GATK_CALLING {

    take:
    sample_cram: Channel<SampleCram>
    reads_grouped: Channel<SampleReads>
    genome: Value<PreparedGenome>
    cohort_size: Value<Integer>
    
    main: 

    // Calculate number of chunks for parallel calling
    included_bases = genome
        .map { prepared ->
            prepared.reference_bases.text.trim().toLong()
        }

    n_chunks = cohort_size
        .combine(included_bases)
        .map { cohort_size_value, included_bases_value ->

            def total_bases =
                cohort_size_value * included_bases_value

            def bases_per_chunk =
                params.hc_target_sample_bases

            Math.max(
                1,
                ((total_bases + bases_per_chunk - 1) / bases_per_chunk) as int
            )
        }

    /* 
        Find and validate any pre-existing GVCFs
    */
    
    // Use existing gvcfs if they are present and the option is set
    if (params.use_existing_gvcf) {

        existing_gvcf = sample_cram
            .map { sample ->

                def vcf = file(
                    "${params.gvcf_store}/${sample.sample}.g.vcf.gz"
                )
                def index = file("${vcf}.tbi")

                record(
                    sample: sample.sample,
                    vcf: vcf,
                    index: index
                )
            }
            .filter { candidate ->
                candidate.vcf.exists() &&
                    candidate.index.exists()
            }

        if (!params.skip_gvcf_validation) {

            gvcf_validation_input = reads_grouped
                .join(
                    existing_gvcf,
                    by: 'sample'
                )

            validation_results = VALIDATE_GVCF(
                gvcf_validation_input,
                genome
            )
                .map { result ->

                    def status = result.status.trim()

                    if (
                        status != 'PASS' &&
                        status != 'FAIL'
                    ) {
                        throw new IllegalStateException(
                            "Unexpected GVCF validation status for " +
                            "${result.sample}: '${status}'"
                        )
                    }

                    record(
                        sample: result.sample,
                        vcf: result.vcf,
                        index: result.index,
                        status: status
                    )
                }

            validated_gvcf = validation_results
                .filter { result ->
                    result.status == 'PASS'
                }
                .map { result ->
                    record(
                        sample: result.sample,
                        vcf: result.vcf,
                        index: result.index
                    )
                }

            validation_results
                .filter { result ->
                    result.status == 'FAIL'
                }
                .map { result ->
                    result.sample
                }
                .unique()
                .collect()
                .subscribe { failed_samples ->
                    if (!failed_samples.isEmpty()) {
                        log.warn(
                            "GVCF validation failed for " +
                                "${failed_samples.size()} sample(s): " +
                                failed_samples.join(', ')
                        )
                    }
                }

        } else {
            validated_gvcf = existing_gvcf
        }

    } else {
        validated_gvcf = channel.empty()
    }

    // Retain crams only for samples without a validated existing GVCF
    gvcf_done = validated_gvcf
        .map { result ->
            result.sample
        }
        .collect()

    cram_for_hc = sample_cram
        .combine(
            done_samples: gvcf_done
        )
        .filter { sample ->
            !sample.done_samples.contains(sample.sample)
        }
        .map { sample ->
            record(
                sample: sample.sample,
                cram: sample.cram,
                crai: sample.crai
            )
        }

    /* 
       Create groups of genomic intervals for parallel calling
       This is done based on CRAI indexes similar to goleft indexsplit
       This is a fast but coarse way of assessign workload
    */

    crai_workload_inputs = sample_cram
        .map { sample ->
            sample.crai
        }
        .collect()

    long_intervals = SUBSET_BED_TO_INTERVALS(
        genome.map { prepared -> prepared.long_bed},
        genome.map { prepared -> prepared.genome_bed},
        genome
    )

    /*
    * For long contigs, project workload only across included calling
    * intervals on long contigs.
    */
    long_projection_genome = genome
        .combine(
            intervals: long_intervals
        )
        .map { prepared ->
            record(
                ref_genome: prepared.ref_genome,
                genome_index_files: prepared.genome_index_files,
                genome_bed: prepared.intervals
            )
        }

    /*
    * For short contigs, project workload across the complete short-contig
    * territory for GenomicsDB compatibility.
    */
    short_projection_genome = genome
        .map { prepared ->
            record(
                ref_genome: prepared.ref_genome,
                genome_index_files: prepared.genome_index_files,
                genome_bed: prepared.short_bed
            )
        }

    long_workload = PROJECT_LONG(
        crai_workload_inputs,
        long_projection_genome
    )

    short_workload = PROJECT_SHORT(
        crai_workload_inputs,
        short_projection_genome
    )

    hc_workload_beds = long_workload
        .combine(short_workload)
        .map { long_bed, short_bed ->
            [long_bed, short_bed]
        }

    combined_workload = COMBINE_WORKLOAD_BEDS(
        hc_workload_beds,
        genome,
        false,
        0,
        4,
        'sum'
    )

    interval_chunks = CREATE_INTERVAL_CHUNKS_HC(
        combined_workload,
        n_chunks,
        params.split_large_intervals
    )

    // Convert the combined workload-balanced chunks into a single
    // genome-wide chunk group.
    hc_chunks = interval_chunks
        .flatMap { chunks ->

            chunks.beds
                .toSorted { bed ->
                    bed.name
                }
                .withIndex()
                .collect { bed, order ->

                    def interval_index = chunks.tbis.find { tbi ->
                        tbi.name == "${bed.name}.tbi"
                    }

                    if (interval_index == null) {
                        throw new IllegalStateException(
                            "Missing index for ${bed.name}"
                        )
                    }

                    def hc_id = bed.name.replaceFirst(/\.bed\.gz$/, '')

                    def jc_number =
                        ((order / params.hc_chunks_per_jc) as Integer) + 1

                    record(
                        jc_id: "genome_jc_${jc_number}",
                        hc_id: hc_id,
                        order: order,
                        interval_bed: bed,
                        interval_index: interval_index
                    )
                }
        }

    // Build one larger JC BED per batch from its HC BEDs.
    jc_bed_inputs = hc_chunks
        .map { chunk ->
            tuple(
                chunk.jc_id,
                chunk.interval_bed
            )
        }
        .groupBy()
        .map { jc_id, hc_files ->
            record(
                jc_id: jc_id,
                hc_files: hc_files
            )
        }

    jc_beds = CREATE_JC_BED_FROM_HC(
        jc_bed_inputs,
        genome
    )

    existing_batch_gvcf = jc_beds
        .map { interval ->
            record(
                jc_id: interval.jc_id
            )
        }
        .combine(validated_gvcf)
        .map { jc, gvcf ->
            record(
                jc_id: jc.jc_id,
                sample: gvcf.sample,
                vcf: gvcf.vcf,
                index: gvcf.index
            )
        }

    /* 
        Single sample calling with HaplotypeCaller
    */

    sample_intervals = hc_chunks
        .combine(cram_for_hc)
        .map { pair ->

            def chunk = pair[0]
            def sample = pair[1]

            record(
                sample: sample.sample,
                jc_id: chunk.jc_id,
                hc_id: chunk.hc_id,
                order: chunk.order,
                interval_bed: chunk.interval_bed,
                interval_index: chunk.interval_index,
                cram: sample.cram,
                crai: sample.crai
            )
        }

    hc_results = HAPLOTYPECALLER(
        sample_intervals,
        genome
    )

    /*
    * Concatenate HC outputs for each sample × JC batch.
    *
    * HAPLOTYPECALLER.out.gvcf_intervals:
    * tuple(sample, jc_id, hc_id, n_hc_in_jc, order, gvcf, tbi)
    */
    jc_concat_plan = hc_results
        .map { result ->
            tuple(
                record(
                    sample: result.sample,
                    jc_id: result.jc_id
                ),
                result
            )
        }
        .groupBy()
        .map { key, files ->

            def ordered = files
                .toSorted { result ->
                    result.order
                }

            record(
                outname: "${key.sample}.${key.jc_id}",
                sample: key.sample,
                jc_id: key.jc_id,
                variants: ordered.collect { result ->
                    result.gvcf
                },
                indices: ordered.collect { result ->
                    result.gvcf_index
                }
            )
        }

    jc_concat_results = CONCAT_GVCFS_JC(
        jc_concat_plan
    )

    jc_concat_ids = jc_concat_plan
        .map { plan ->
            record(
                outname: plan.outname,
                sample: plan.sample,
                jc_id: plan.jc_id
            )
        }

    new_batch_gvcf = jc_concat_results
        .join(
            jc_concat_ids,
            by: 'outname'
        )
        .map { result ->
            record(
                jc_id: result.jc_id,
                sample: result.sample,
                vcf: result.vcf,
                index: result.index
            )
        }

    /*
    * Release a JC batch for processing once its expected number of samples is ready.
    * Keep each sample, gVCF and index together while grouping.
    */
    gvcf_intervals = new_batch_gvcf
        .mix(existing_batch_gvcf)
        .map { gvcf ->
            tuple(
                gvcf.jc_id,
                record(
                    sample: gvcf.sample,
                    vcf: gvcf.vcf,
                    index: gvcf.index
                )
            )
        }
        .groupBy()
        .map { jc_id, files ->

            def ordered = files
                .toSorted { result ->
                    result.sample
                }

            def samples = ordered.collect { result ->
                result.sample
            }

            if (samples.size() != samples.toSet().size()) {
                throw new IllegalStateException(
                    "Duplicate sample in JC batch ${jc_id}: ${samples}"
                )
            }

            record(
                jc_id: jc_id,
                gvcfs: ordered.collect { result ->
                    result.vcf
                },
                indices: ordered.collect { result ->
                    result.index
                }
            )
        }
        .join(
            jc_beds,
            by: 'jc_id'
        )
        .map { result ->
            record(
                jc_id: result.jc_id,
                interval_bed: result.interval_bed,
                interval_index: result.interval_index,
                gvcfs: result.gvcfs,
                indices: result.indices
            )
        }

    genomicsdb_results = GENOMICSDB_IMPORT(
        gvcf_intervals,
        genome,
        cohort_size
    )

    joint_genotype_results = JOINT_GENOTYPE(
        genomicsdb_results,
        genome,
        cohort_size
    )

    // Normalise back to cannonical variant type
    joint_variants = joint_genotype_results
        .map { result ->
            record(
                interval_id: result.jc_id,
                interval_bed: result.interval_bed,
                interval_index: result.bed_tbi,
                variants: result.bcf,
                variants_index: result.bcf_index
            )
        }

    
    /* 
    * Extra publishing
    */
    
    // Per-sample gVCFs
    if (params.output_gvcf) {

        whole_sample_concat = hc_results
            .map { result ->
                tuple(
                    result.sample,
                    record(
                        hc_id: result.hc_id,
                        order: result.order,
                        gvcf: result.gvcf,
                        gvcf_index: result.gvcf_index
                    )
                )
            }
            .groupBy()
            .map { sample, pieces ->

                def ordered = pieces
                    .toSorted { piece ->
                        piece.order
                    }

                record(
                    outname: sample,
                    variants: ordered.collect { piece ->
                        piece.gvcf
                    },
                    indices: ordered.collect { piece ->
                        piece.gvcf_index
                    }
                )
            }

        new_gvcf = CONCAT_GVCFS_SAMPLE(
            whole_sample_concat
        )

    } else {

        new_gvcf = channel.empty()
    }


    // Unfiltered cohort BCF
    merged_unfiltered_variants = channel.empty()

    if (params.output_unfiltered_bcf) {

        variants_to_merge = joint_variants
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


    emit:
    new_gvcf: Channel<VariantFile> = new_gvcf
    variants: Channel<VariantChunk> = joint_variants
    merged_unfiltered_variants: Channel<VariantFile> = merged_unfiltered_variants
}