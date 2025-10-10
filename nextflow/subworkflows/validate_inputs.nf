/*
    Validate inputs
*/

//// import modules
include { VALIDATE_FASTQ                        } from '../modules/validate_fastq'
include { REPAIR_FASTQ                          } from '../modules/repair_fastq'

workflow VALIDATE_INPUTS {

    take:
    ch_sample_names
    ch_reads
    ch_existing_cram
    ch_existing_gvcf

    main: 

    /* 
        Validate fastq and extract read groups
    */
    VALIDATE_FASTQ (
        ch_reads
    )

    // Convert stdout to a string for status (PASS or FAIL)
    VALIDATE_FASTQ.out.fastq_with_status
        .map { sample, lib, read1, read2, stdout -> [ sample, lib, read1, read2, stdout.trim() ] }
        .branch { sample, lib, read1, read2, status ->
            fail: status == 'FAIL'
            pass: status == 'PASS'
        }
        .set { validation_routes }

    // Print a warning if any samples fail validation and need to be repaired
    validation_routes.fail
    .map { sample, lib, read1, read2, _ -> sample } 
    .unique()
    .collect()
    .map { fails ->
        if (fails && fails.size() > 0)
        log.warn "Repairing malformed FASTQs for ${fails.size()} sample(s): ${fails.join(', ')}"
        true
    }
    .set { _warn_done }  // force evaluation

    // Repair any fastqs that failed validation 
    REPAIR_FASTQ(
        validation_routes.fail.map { sample, lib, read1, read2, _ -> [sample, lib, read1, read2] }
    )

    // Join repaired fastqs back into validated fastqs
    validation_routes.pass.map { sample, lib, read1, read2, _ -> [sample, lib, read1, read2] }
        .mix( REPAIR_FASTQ.out.fastq )
        .set { ch_validated_fastq }

    // Validate existing CRAMS
    // TODO:: to pass validation the .cram and .crai must exist AND the readgroups must contain all FASTQ readgroups for that sample
    ch_existing_cram
        .set{ ch_validated_cram }

    // Validate GVCFs
    // TODO:: to pass validation the .gvcf and tbi must exist AND the comment line must contain all FASTQ readgroups for that sample
    ch_existing_gvcf
        .set{ ch_validated_gvcf }
    
    emit: 
    validated_fastq = ch_validated_fastq
    validated_cram = ch_validated_cram
    validated_gvcf = ch_validated_gvcf

}

