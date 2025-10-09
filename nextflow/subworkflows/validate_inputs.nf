/*
    Validate inputs
*/

//// import modules
include { VALIDATE_FASTQ                        } from '../modules/validate_fastq'
include { REPAIR_FASTQ                          } from '../modules/repair_fastq'

workflow VALIDATE_INPUTS {

    take:
    ch_reads
    ch_existing_cram
    ch_existing_gvcf

    main: 

    // TODO: Validate GVCFs

    // TODO: Validate CRAMS first

    // keys for reads
    ch_reads
        .map { s, lib, r1, r2 -> tuple(s, [s, lib, r1, r2]) }
        .set { ch_reads_by_sample }

    // keys for done samples (only if both files exist)
    ch_existing_cram
        .map    { s, cram, crai -> tuple(s, true) }
        .set { ch_done_keys }

    // left anti-join: keep reads with no match in done keys
    ch_reads_by_sample
        .join(ch_done_keys, remainder: true)
        .filter { sample, payload, doneFlag -> doneFlag == null }
        .map    { sample, payload, _ -> payload }
        .set    { ch_reads_to_validate }

    /* 
        Validate fastq and extract read headers
    */
    VALIDATE_FASTQ (
        ch_reads_to_validate
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
        .set { ch_all_fixed_fastq }
    
    emit: 
    reads_to_map = ch_all_fixed_fastq
    validated_cram = ch_existing_cram

}

