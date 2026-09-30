/*
    Genotype samples using GATK
*/

//// import modules
include { VALIDATE_GVCF                                          } from '../modules/validate_gvcf/validate_gvcf'
include { COMBINE_MOSDEPTH_EVENTS                                } from '../modules/combine_mosdepth_events/combine_mosdepth_events'
include { SPLIT_BED_BY_CHR                                       } from '../modules/split_bed_by_chr/split_bed_by_chr' 
include { CREATE_JC_BED_FROM_HC                                  } from '../modules/create_jc_bed_from_hc/create_jc_bed_from_hc' 
include { CREATE_INTERVAL_CHUNKS as CREATE_IC_LONG               } from '../modules/create_interval_chunks/create_interval_chunks'
include { CREATE_INTERVAL_CHUNKS as CREATE_IC_SHORT              } from '../modules/create_interval_chunks/create_interval_chunks'
include { HAPLOTYPECALLER                                        } from '../modules/haplotypecaller/haplotypecaller'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import/genomicsdb_import' 
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype/joint_genotype' 
include { STAGE_GVCF                                             } from '../modules/stage_gvcf/stage_gvcf'
include { CONCAT_VCFS as CONCAT_GVCFS_JC                         } from '../modules/concat_vcfs/concat_vcfs'
include { CONCAT_VCFS as CONCAT_GVCFS_SAMPLE                     } from '../modules/concat_vcfs/concat_vcfs'
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                  } from '../modules/concat_vcfs/concat_vcfs' 

workflow GATK_CALLING {

    take:
    ch_sample_names
    ch_sample_cram
    ch_reads_grouped
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_genotype
    ch_genome_depths
    ch_long_bed
    ch_short_bed
    
    main: 

    /* 
        Find and validate any pre-existing GVCFs
    */
    
    // Use existing gvcfs if they are present and the option is set
    if( params.use_existing_gvcf ) {
        ch_sample_names
            .map { sample ->
                def gvcf = file("${params.gvcf_store}/${sample}.g.vcf.gz")
                def tbi = file("${gvcf}.tbi")
                tuple(sample, gvcf, tbi)
            }
            .filter { _sample, gvcf, tbi -> gvcf.exists() && tbi.exists() }
            .set { ch_existing_gvcf }


        // Validate gvcf files by default
        if( !params.skip_gvcf_validation ) {

            ch_reads_grouped
                .join(ch_existing_gvcf, by: 0)
                .set { ch_gvcf_validation_input }

            VALIDATE_GVCF (
                ch_gvcf_validation_input,
                ch_genome_indexed
            )

            // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
            VALIDATE_GVCF.out.status
                .map { sample, stdout -> [ sample, stdout.trim() ] }
                .join( ch_existing_gvcf, by: 0 )
                .map { sample, status, gvcf, tbi -> [ sample, gvcf, tbi, status ] }
                .branch {  _sample, _gvcf, _tbi, status ->
                    fail: status == 'FAIL'
                    pass: status == 'PASS'
                    invalid: true
                }
                .set { gvcf_validation_routes }

            // Fail loudly if validator emits anythign other that pass/fail
            gvcf_validation_routes.invalid
                .map { sample, _gvcf, _tbi, status ->
                    throw new IllegalStateException(
                        "Unexpected GVCF validation status for " +
                        "${sample}: '${status}'"
                    )
                }
                .set { _invalid_gvcf_status }

            // Compatible existing GVCFs.
            gvcf_validation_routes.pass
                .map { sample, gvcf, tbi, _status -> [ sample, gvcf, tbi ] } 
                .set { ch_validated_gvcf }
                
            // Print warning if any gvcf files exist but fail validation
            gvcf_validation_routes.fail
                .map {  sample, _gvcf, _tbi, _status -> sample } 
                .unique()
                .collect()
                .subscribe { failed_samples ->
                    if (failed_samples) {
                        log.warn(
                            "GVCF validation failed for " +
                            "${failed_samples.size()} sample(s): " +
                            failed_samples.join(', ')
                        )
                    }
                }

        } else {
          // Skip validation, assume all existing gvcfs are compatible
          ch_validated_gvcf = ch_existing_gvcf 
        }

        // Set of samples that do not need single-sample variant calling.
        ch_validated_gvcf
            .map { sample, _gvcf, _tbi -> sample}
            .toList()
            .map { ids -> ids as Set } 
            .set { ch_gvcf_done }

    } else{
        ch_gvcf_done = channel.value([] as Set)
        ch_validated_gvcf = channel.empty()
    }

    // Retain crams only for samples without a validated existing GVCF
    ch_sample_cram
        .combine(ch_gvcf_done)  
        .filter { sample, _gvcf, _tbi, doneSet -> !(doneSet as Set).contains(sample) }
        .map {  sample, gvcf, tbi, _doneSet -> tuple( sample, gvcf, tbi) }
        .set { ch_cram_for_hc }


    /* 
        Create interval chunks
    */
    
    // First split bed by chr    
    // These must already contain full-contig coordinates.
    SPLIT_BED_BY_CHR(ch_long_bed.first())

    // Create interval chunks for long bed  (keep_whole_contigs = false)
    CREATE_IC_LONG(
        SPLIT_BED_BY_CHR.out.per_chr_beds.flatten(),
        COMBINE_MOSDEPTH_EVENTS.out.rle.first(),
        params.hc_bases_per_chunk,
        params.min_interval_gap,
        false
    )

    // Create interval chunks for short bed (keep_whole_contigs = true)
    CREATE_IC_SHORT(
        ch_short_bed.first(),
        COMBINE_MOSDEPTH_EVENTS.out.rle.first(),
        params.hc_bases_per_chunk,
        0,
        true
    )

    // Mix long and short contigs back together
    CREATE_IC_LONG.out.interval_bed
        .map { include_bed, beds, tbis -> tuple("long_${include_bed}", beds, tbis)}
        .mix( CREATE_IC_SHORT.out.interval_bed.map { _include_bed, beds, tbis -> tuple("short", beds, tbis) })
        .map { region_id, beds, tbis ->
            def bedList = beds instanceof List ? beds : [beds]
            def tbiByName = (tbis instanceof List ? tbis : [tbis])
                .collectEntries { tbi -> [(tbi.name): tbi] }

            def chunks = bedList.collect { bed ->
                def tbi = tbiByName["${bed.name}.tbi"]

                assert tbi != null :
                    "Missing index for ${region_id}: ${bed.name}"

                def hc_id = "${region_id}_${bed.name.replaceFirst(/\.bed\.gz$/, '')}"
                tuple(region_id, hc_id, bed, tbi)
            }.sort { left, right -> left[1] <=> right[1] }

            tuple(region_id, chunks)
        }
        .set { ch_hc_chunks_by_region }

    // Build one larger JC BED per batch from its HC BEDs.

    // Plan JC batches independently for each long contig or
    // for the short-contig group. No genome-wide toList().
    ch_hc_chunks_by_region
        .flatMap { region_id, chunks ->
            chunks
                .collate(params.hc_chunks_per_jc)
                .withIndex()
                .collect { group, index ->
                    def jc_id =
                        "${region_id}_jc_${String.format('%05d', index + 1)}"

                    tuple(jc_id, region_id, group)
                }
        }
        .set { ch_jc_plan }

    CREATE_JC_BED_FROM_HC(
        ch_jc_plan.map { jc_id, _region_id, group ->
            tuple(jc_id, group.collect { _region, _hc_id, bed, _tbi -> bed })
        },
        ch_genome_indexed
    )


    /* 
        Single sample calling with HaplotypeCaller
    */

    // ch_jc_plan: tuple(jc_id, region_id, group)
    // group entries: tuple(region_id, hc_id, bed, tbi)
    ch_jc_plan
        .flatMap { jc_id, _region_id, group ->
            group.withIndex().collect { chunk, order ->
                tuple(jc_id, chunk[1], group.size(), order, chunk[2], chunk[3])
            }
        }
        .combine(ch_cram_for_hc)
        .map { jc_id, hc_id, n_hc_in_jc, order, bed, tbi,
            sample, cram, crai ->
            tuple(
                sample, jc_id, hc_id, n_hc_in_jc, order,
                bed, tbi, cram, crai
            )
        }
        .set { ch_sample_intervals }

    HAPLOTYPECALLER(
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_genotype
    )

    /*
    * Concatenate HC outputs for each sample × JC batch.
    *
    * HAPLOTYPECALLER.out.gvcf_intervals:
    * tuple(sample, jc_id, hc_id, n_hc_in_jc, order, gvcf, tbi)
    */
    HAPLOTYPECALLER.out.gvcf_intervals
        .map { sample, jc_id, _hc_id, n_hc_in_jc, _order, gvcf, tbi ->
            tuple(
                groupKey([sample, jc_id], n_hc_in_jc),
                tuple(gvcf, tbi)
            )
        }
        .groupTuple()
        .map { key, files ->
            def (sample, jc_id) = key.getGroupTarget()

            tuple(
                "${sample}.${jc_id}",
                sample,
                jc_id,
                files.collect { file_pair -> file_pair[0] },
                files.collect { file_pair -> file_pair[1] }
            )
        }
        .set { ch_jc_concat_plan }

    // Keep identifiers for restoring sample and jc_id after concat.
    ch_jc_concat_plan
        .map { outname, sample, jc_id, _gvcfs, _tbis ->
            tuple(outname, sample, jc_id)
        }
        .set { ch_jc_concat_ids }

    // Existing CONCAT_VCFS input: tuple(outname, gvcfs, tbis).
    ch_jc_concat_plan
        .map { outname, _sample, _jc_id, gvcfs, tbis ->
            tuple(outname, gvcfs, tbis)
        }
        .set { ch_jc_concat_input }

    CONCAT_GVCFS_JC(ch_jc_concat_input)

    // tuple(jc_id, sample, batch_gvcf, batch_tbi)
    CONCAT_GVCFS_JC.out.vcf
        .join(ch_jc_concat_ids, by: 0)
        .map { _outname, gvcf, tbi, sample, jc_id ->
            tuple(jc_id, sample, gvcf, tbi)
        }
        .set { ch_new_batch_gvcf }

    /*
    * Add each validated whole-sample gVCF to every JC batch.
    * The matching JC BED will restrict the import territory.
    */
    CREATE_JC_BED_FROM_HC.out.interval_bed
        .map { jc_id, _bed, _bed_tbi -> jc_id }
        .combine(ch_validated_gvcf)
        .map { jc_id, sample, gvcf, tbi ->
            tuple(jc_id, sample, gvcf, tbi)
        }
        .set { ch_existing_batch_gvcf }

    // Reuse this value channel both for grouping and resource scaling.
    ch_cohort_size = ch_sample_names.unique().count()

    /*
    * Release a JC batch for processing once its expected number of samples is ready.
    * Keep each sample, gVCF and index together while grouping.
    */
    ch_new_batch_gvcf
        .mix(ch_existing_batch_gvcf)
        .combine(ch_cohort_size)
        .map { jc_id, sample, gvcf, tbi, n_samples ->
            tuple(
                groupKey(jc_id, n_samples),
                tuple(sample, gvcf, tbi)
            )
        }
        .groupTuple()
        .map { key, files ->
            def jc_id = key.getGroupTarget()
            def ordered = files.sort { a, b -> a[0] <=> b[0] }
            def samples = ordered.collect { entry -> entry[0] }

            assert samples.size() == samples.toSet().size() :
                "Duplicate sample in JC batch ${jc_id}: ${samples}"

            tuple(
                jc_id,
                ordered.collect { entry -> entry[1] },
                ordered.collect { entry -> entry[2] }
            )
        }
        .join(CREATE_JC_BED_FROM_HC.out.interval_bed, by: 0)
        .map { jc_id, gvcfs, gvcf_tbis, jc_bed, jc_bed_tbi ->
            // Matches the existing GENOMICSDB_IMPORT five-field input.
            tuple(jc_id, jc_bed, jc_bed_tbi, gvcfs, gvcf_tbis)
        }
        .set { ch_gvcf_interval }

    GENOMICSDB_IMPORT(
        ch_gvcf_interval,
        ch_genome_indexed,
        ch_cohort_size
    )

    JOINT_GENOTYPE(
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_genotype,
        ch_cohort_size
    )

    
    /* 
    * Extra publishing
    */
    
    // Per-sample gvcfs
    ch_new_gvcf = channel.empty()
    if ( params.output_gvcf ){
        HAPLOTYPECALLER.out.gvcf_intervals
            .map { sample, _jc_id, hc_id, _n_hc_in_jc, _order, gvcf, tbi ->
                tuple(sample, tuple(hc_id, gvcf, tbi))
            }
            .groupTuple(by: 0)
            .map { sample, pieces ->
                // Keep each gVCF paired with its index.
                // Do not assume HC completion order is genomic order.
                tuple(
                    sample,
                    pieces.collect { piece -> piece[1] },
                    pieces.collect { piece -> piece[2] }
                )
            }
            .set { ch_whole_sample_to_concat }

        CONCAT_GVCFS_SAMPLE(ch_whole_sample_to_concat)

        CONCAT_GVCFS_SAMPLE.out.vcf
            .set { ch_new_gvcf }

        ch_validated_gvcf
            .mix(ch_new_gvcf)
            .set { ch_sample_gvcf }

        STAGE_GVCF(ch_sample_gvcf)
    }

    // unfiltered cohort vcf
    ch_merged_unfiltered_vcf = channel.empty()
    if ( params.output_unfiltered_vcf ){

        // TODO: Make this output seperate files for each variant type
        JOINT_GENOTYPE.out.vcf
            .map { _interval_chunk, _interval_bed, _bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .map { _type, vcf, tbi -> tuple('all', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        CONCAT_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )

        CONCAT_UNFILTERED_VCFS.out.vcf
            .set { ch_merged_unfiltered_vcf }
    }

    emit: 
    new_gvcf = ch_new_gvcf
    vcf = JOINT_GENOTYPE.out.vcf
    merged_unfiltered_vcf = ch_merged_unfiltered_vcf
}