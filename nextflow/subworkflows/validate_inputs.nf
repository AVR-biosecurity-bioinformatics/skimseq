/*
    Validate inputs
*/

//// import modules
include { VALIDATE_FASTQ                        } from '../modules/validate_fastq'
include { REPAIR_FASTQ                          } from '../modules/repair_fastq'

workflow VALIDATE_INPUTS {

    take:
    ch_reads

    main: 

    /* 
        Validate fastq and extract read headers
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
        .set { ch_all_fixed_fastq }
    
    emit: 
    reads_to_map = ch_all_fixed_fastq
}

