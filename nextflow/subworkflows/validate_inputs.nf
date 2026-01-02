/*
    Validate inputs
*/

//// import modules
include { VALIDATE_FASTQ                       } from '../modules/validate_fastq'
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

    emit: 
    validated_fastq = ch_validated_fastq
    rg_to_validate = ch_rg_to_validate

}