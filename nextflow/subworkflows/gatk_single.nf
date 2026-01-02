/*
    Genotype samples using GATK
*/

//// import modules
include { VALIDATE_GVCF                                          } from '../modules/validate_gvcf'
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { COUNT_VCF_DEPTH                                        } from '../modules/count_vcf_depth'
include { CREATE_INTERVAL_CHUNKS_HC                              } from '../modules/create_interval_chunks_hc'
include { PROFILE_HC                                             } from '../modules/profile_hc'
include { STAGE_GVCF                                             } from '../modules/stage_gvcf'

workflow GATK_SINGLE {

    take:
    ch_sample_names
    ch_sample_cram
    ch_rg_to_validate
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed
    ch_dummy_file

    main: 

    /* 
        Find and validate any pre-existing GVCFs
    */
    
    ch_sample_names
        .map { sample ->
            def gvcf = file("output/results/vcf/gvcf/${sample}.g.vcf.gz")
            def tbi = file("${gvcf}.tbi")
            tuple(sample, gvcf, tbi)
        }
        .filter { sample, gvcf, tbi -> gvcf.exists() && tbi.exists() }
        .set { ch_existing_gvcf }

    // Validate GVCFs
    VALIDATE_GVCF (
        ch_rg_to_validate.join(ch_existing_gvcf, by: 0),
        ch_genome_indexed
    )

    // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
    VALIDATE_GVCF.out.status
        .map { sample, stdout -> [ sample, stdout.trim() ] }
        .join( ch_existing_gvcf, by: 0 )
        .map { sample, status, gvcf, tbi -> [ sample, gvcf, tbi, status ] }
        .branch {  sample, gvcf, tbi, status ->
            fail: status == 'FAIL'
            pass: status == 'PASS'
        }
        .set { gvcf_validation_routes }

    gvcf_validation_routes.pass
        .map { sample, gvcf, tbi, status -> [ sample, gvcf, tbi ] } 
        .set { ch_validated_gvcf }
        
    // Print warning if any gvcf files exist but fail validation
    gvcf_validation_routes.fail
        .map {  sample, gvcf, tbi, status -> sample } 
        .unique()
        .collect()
        .map { fails ->
            if (fails && fails.size() > 0)
            log.warn "GVCF file failed validation for ${fails.size()} samples(s): ${fails.join(', ')}"
            true
        }
        .set { _warn_gvcf_done }  // force evaluation

    // Subset the crams to just those that dont already have a GVCF for single sample calling
    ch_validated_gvcf
        .map { sample, gvcf, tbi -> sample}
        .toList()
        .map { ids -> ids as Set } 
        .set { ch_gvcf_done }

    ch_sample_cram
        .combine(ch_gvcf_done)  
        .filter { sample, gvcf, tbi, doneSet -> !(doneSet as Set).contains(sample) }
        .map {  sample, gvcf, tbi, doneSet -> tuple( sample, gvcf, tbi) }
        .set { ch_cram_for_hc }

    /* 
       Create groups of genomic intervals for parallel haplotypecaller
    */

    // Create haplotypecaller intervals on per sample basis
    CREATE_INTERVAL_CHUNKS_HC (
        ch_cram_for_hc,
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        params.hc_bases_per_chunk,
        params.split_large_intervals,
        params.hc_rmdup,
        params.hc_minbq,
        params.hc_minmq,
        params.min_interval_gap
    )
   
    // CREATE_INTERVAL_CHUNKS_HC.out.interval_bed emits: tuple(sample, bed)
    // where `bed` is either a List<Path> or a single Path, so has to be normalised to list
    CREATE_INTERVAL_CHUNKS_HC.out.interval_bed
    .flatMap { sample, beds ->
        // normalize to a list for cases where there are only 1 bed output for a sample
        def lst = (beds instanceof List) ? beds : [ beds ]
        // emit one tuple per bed file
        lst.collect { bed ->
        bed  = bed as Path
        def base = bed.baseName
        def interval_chunk = base.startsWith('_') ? base.substring(1) : base
        tuple(sample, interval_chunk, bed)
        }
    }
    .set { ch_interval_bed_hc }

    // Combine intervals with cram files for genotyping
    ch_interval_bed_hc 
	.combine( ch_cram_for_hc, by: [0, 0] )
    .set { ch_sample_intervals }

    /* 
       Call variants per sample
    */

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_gatk
    )

    if( params.profile_gatk ) {

        // Join back onto cram and gvcf based on first 3 columns
        CALL_VARIANTS.out.log
            .join( ch_sample_intervals.map { sample, interval_chunk, interval_bed,cram, crai -> tuple(sample, interval_chunk, cram, crai) }, by:[0,1] )
            .join( CALL_VARIANTS.out.gvcf_intervals, by:[0,1] )
            .map { sample, interval_chunk, logfile, assembly_regions, cram, crai, gvcf, tbi -> tuple(sample, interval_chunk, cram, crai, gvcf, tbi, logfile, assembly_regions ) }
            .set { ch_for_profile }

        // Profile HC runtimes per Sample x Interval
        PROFILE_HC (
            ch_for_profile,
            ch_genome_indexed
        )

        // Merge and output HC profiles
        PROFILE_HC.out.summary
            .collectFile(
                name: 'hc_profiles.tsv',
                storeDir: "${launchDir}/output/gatk_profiles",
                skip: 1,
                keepHeader: true,
                newLine: false,
                sort: true
            )
    }

    // Merge interval GVCFs by sample
    CALL_VARIANTS.out.gvcf_intervals
        .groupTuple ( by: 0 )
        .set { ch_gvcf_to_merge }

    MERGE_GVCFS (
        ch_gvcf_to_merge.map { sample, interval_chunk, gvcf, tbi -> [ sample, gvcf, tbi ] }
    )

    // Update the newly created gvcf path to the canonical publishdir path to ensure that resume works for further steps
    //MERGE_GVCFS.out.vcf
    //    .map { sample, gvcf, tbi ->
    //        def realgvcf = file("output/results/vcf/gvcf/${sample}.g.vcf.gz")
    //        def realtbi = file("output/results/vcf/gvcf/${sample}.g.vcf.gz.tbi")
    //        tuple(sample, realgvcf, realtbi)
    //    }
    //    .set { ch_new_gvcf_canonical }

    // combine validated existing GVCs with newly created GVCFs for joint calling
    ch_validated_gvcf
      .mix( MERGE_GVCFS.out.vcf )
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