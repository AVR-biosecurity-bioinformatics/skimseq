/*
    Validate inputs
*/

//// import modules
include { VALIDATE_FASTQ                       } from '../modules/validate_fastq'
include { VALIDATE_CRAM                        } from '../modules/validate_cram'
include { VALIDATE_GVCF                        } from '../modules/validate_gvcf'
include { REPAIR_FASTQ                         } from '../modules/repair_fastq'

workflow VALIDATE_INPUTS {

    take:
    ch_sample_names
    ch_reads
    ch_genome_indexed

    main: 

    /* 
        Validate fastq and extract read groups
    */
    VALIDATE_FASTQ (
        ch_reads
    )

    // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
    VALIDATE_FASTQ.out.status
        .map { sample, lib, stdout ->
            def (fcid, lane, platform, status) = stdout.trim().tokenize('\t')
            tuple(sample, lib, fcid, lane, platform, status)
        }
        .join( ch_reads, by:[0,1] )
        .map { sample, lib, fcid, lane, platform, status, read1, read2 -> [ sample, lib, fcid, lane, platform, read1, read2, status ] }
        .branch { sample, lib, fcid, lane, platform, read1, read2, status ->
            fail: status == 'FAIL'
            pass: status == 'PASS'
        }
        .set { fastq_validation_routes }

    // Print a warning if any samples fail validation and need to be repaired
    fastq_validation_routes.fail
    .map { sample, lib, fcid, lane, platform, read1, read2, status -> lib } 
    .unique()
    .collect()
    .map { fails ->
        if (fails && fails.size() > 0)
        log.warn "Repairing malformed FASTQs for ${fails.size()} libraries(s): ${fails.join(', ')}"
        true
    }
    .set { _warn_done }  // force evaluation

    // Repair any fastqs that failed validation 
    REPAIR_FASTQ(
        fastq_validation_routes.fail.map { sample, lib, fcid, lane, platform, read1, read2, status -> [sample, lib, fcid, lane, platform, read1, read2] }
    )

    // Join repaired fastqs back into validated fastqs
    fastq_validation_routes.pass.map { sample, lib, fcid, lane, platform, read1, read2, status -> [sample, lib, fcid, lane, platform, read1, read2] }
        .mix( REPAIR_FASTQ.out.fastq )
        .set { ch_validated_fastq }

    // Create a rg to validate channel to be used for cram and gvcf validateion
    ch_validated_fastq
        //.map { sample, lib, fcid, lane, platform, read1, read2 -> [sample, lib, fcid, lane, platform] }
        .groupTuple(by: 0)
        .map {sample, lib, fcid, lane, platform, read1, read2 ->
                // nested list of RG metadata per library
                def rg_list = (0..<lib.size()).collect { i ->
                    [ sample, lib[i], fcid[i], lane[i], platform[i] ]
                }

                // nested list of read1 fastqs (one per library)
                def r1_list = (0..<lib.size()).collect { i ->
                    read1[i]
                }

                // nested list of read2 fastqs (one per library)
                def r2_list = (0..<lib.size()).collect { i ->
                    read2[i]
                }

                // emit sample with all three nested lists
                tuple(sample, rg_list, r1_list, r2_list)
        }
        .set { ch_rg_to_validate }

    /* 
        Find and validate any pre-existing crams
        To pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    */

         ch_sample_names
        .map { sample ->
            def cram = file("output/results/cram/${sample}.cram")
            def crai = file("${cram}.crai")
            tuple(sample, cram, crai)
        }
        .filter { sample, cram, crai -> cram.exists() && crai.exists() }
        .set { ch_existing_cram }

    // Validate cram files
    VALIDATE_CRAM (
        ch_rg_to_validate.join(ch_existing_cram, by: 0),
        ch_genome_indexed
    )

    // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
    VALIDATE_CRAM.out.status
        .map { sample, stdout -> [ sample, stdout.trim() ] }
        .join( ch_existing_cram, by: 0 )
        .map { sample, status, cram, crai -> [ sample, cram, crai, status ] }
        .branch {  sample, cram, crai, status ->
            fail: status == 'FAIL'
            pass: status == 'PASS'
        }
        .set { cram_validation_routes }

    cram_validation_routes.pass
        .map { sample, cram, crai, status -> [ sample, cram, crai ] } 
        .set { ch_validated_cram }
        
    // Print warning if any cram files exist but fail validation
    cram_validation_routes.fail
        .map {  sample, cram, crai, status -> sample } 
        .unique()
        .collect()
        .map { fails ->
            if (fails && fails.size() > 0)
            log.warn "CRAM file failed validation for ${fails.size()} samples(s): ${fails.join(', ')}"
            true
        }
        .set { _warn_cram_done }  // force evaluation


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

    // TODO:: to pass validation the the comment line must contain all FASTQ readgroups for that sample
    //ch_existing_gvcf
    //    .set{ ch_validated_gvcf }
    
    emit: 
    validated_fastq = ch_validated_fastq
    validated_cram = ch_validated_cram
    validated_gvcf = ch_validated_gvcf

}