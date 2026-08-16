/*
    Genotype samples using GATK
*/

//// import modules
include { VALIDATE_GVCF                                          } from '../modules/validate_gvcf/validate_gvcf'
include { HAPLOTYPECALLER                                        } from '../modules/haplotypecaller/haplotypecaller'
include { CONCAT_VCFS as CONCAT_GVCFS                            } from '../modules/concat_vcfs/concat_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_HC    } from '../modules/create_interval_chunks/create_interval_chunks'
include { STAGE_GVCF                                             } from '../modules/stage_gvcf/stage_gvcf'

workflow GATK_SINGLE {

    take:
    ch_sample_names
    ch_sample_cram
    ch_reads_grouped
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_genotype
    ch_long_bed
    ch_short_bed
    ch_read_counts
    
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

    ch_sample_cram
        .combine(ch_gvcf_done)  
        .filter { sample, _gvcf, _tbi, doneSet -> !(doneSet as Set).contains(sample) }
        .map {  sample, gvcf, tbi, _doneSet -> tuple( sample, gvcf, tbi) }
        .set { ch_cram_for_hc }

    /* 
       Create groups of genomic intervals for parallel haplotypecaller
    */

    // Create haplotypecaller intervals on per sample basis
    CREATE_INTERVAL_CHUNKS_HC (
        ch_read_counts,
        ch_genome_indexed,
        ch_include_bed.first(),
        params.hc_bases_per_chunk,
        params.min_interval_gap,
        params.split_large_intervals,
        "false"
    )
   
    // CREATE_INTERVAL_CHUNKS_HC.out.interval_bed emits: tuple(sample, bed)
    // where `bed` is either a List<Path> or a single Path, so has to be normalised to list
    CREATE_INTERVAL_CHUNKS_HC.out.interval_bed
        .flatMap { sample, beds, tbis  ->
            // normalize to a list for cases where there are only 1 bed output for a sample
            def bedList = (beds instanceof List) ? beds : [beds]
            def tbiList = (tbis instanceof List) ? tbis : [tbis]

            assert bedList.size() == tbiList.size() :
            "Mismatch for ${sample}: beds=${bedList.size()} tbis=${tbiList.size()}"

            // Count number of intervals expected for that sample
            def n_intervals = bedList.size()

            // emit one tuple per bed file
            (0..<bedList.size()).collect { i ->
                def bed = bedList[i] as Path
                def tbiPath = tbiList[i]
                def base = bed.getFileName().toString()
                base = base.replaceFirst(/\.gz$/, '')
                base = base.replaceFirst(/\.bed$/, '')
                def interval_hash = base.startsWith('_') ? base.substring(1) : base
                tuple(sample, interval_hash, n_intervals, bed, tbiPath)
            }
        }
        .set { ch_interval_bed_hc }

    // Combine intervals with cram files for genotyping
    ch_interval_bed_hc
        .combine(ch_cram_for_hc, by: 0)
        .map { sample, interval_hash, n_intervals, bed, tbi, cram, crai ->
            tuple(sample, interval_hash, n_intervals, bed, tbi, cram, crai)
        }
        .set { ch_sample_intervals }

    // TODO: are number of intervals per sample the same?
    // Can i just count how many intervals there are per sample in ch_sample_intervals? - or does this require waiting too
    // Or output nchunks from CREATE_INTERVAL_CHUNKS

    /* 
       Call variants per sample
    */

    // call variants for single samples across intervals
    HAPLOTYPECALLER (
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_genotype
    )

    // Grouping by sample, nchunks allows early per-sample merge rather than waiting for all HAPLOTYPECALLER to finish
    HAPLOTYPECALLER.out.gvcf_intervals
        // Add expected group size so each sample emits once all interval GVCFs are complete
        .map { sample, _interval_hash, n_intervals, gvcf, tbi ->
            tuple(groupKey(sample, n_intervals), gvcf, tbi)
        }
        // Group interval GVCFs by sample, emitting early when n_intervals have arrived
        .groupTuple()
        // Emit sample with grouped GVCFs and indexes for concatenation
        .map { key, gvcfs, tbis ->
            tuple(
                key.getGroupTarget(),
                gvcfs,
                tbis
            )
        }
        // channel for per-sample GVCF concatenation
        .set { ch_gvcf_to_concat }


    CONCAT_GVCFS (
        ch_gvcf_to_concat
    )

    // combine validated existing GVCs with newly created GVCFs for joint calling
    ch_validated_gvcf
      .mix( CONCAT_GVCFS.out.vcf )
      .distinct { it[0] }      // dedupe by sample if needed
      .set{ ch_sample_gvcf }

    // Helper process to publish to output directory. 
    // NOTE: This process (using deep caching) is necessary to avoid violating cache of later steps when inputs switch to existing gvcf on resume
    STAGE_GVCF(
        ch_sample_gvcf
    )

    emit: 
    gvcf = STAGE_GVCF.out.gvcf
}