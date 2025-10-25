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

    main: 

    /* 
        Validate fastq and extract read groups
    */
    VALIDATE_FASTQ (
        ch_reads
    )

    // Convert stdout to a string for status (PASS or FAIL), and join to initial reads
    VALIDATE_FASTQ.out.status
        .map { sample, lib, stdout -> [ sample, lib, stdout.trim() ] }
        .join( ch_reads, by:[0,1] )
        .map { sample, lib, status, read1, read2 -> [ sample, lib, read1, read2, status ] }
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

    /* 
        Find and validate any pre-existing crams
    */

    ch_sample_names
        .map { s ->
            def cram = file("output/results/cram/${s}.cram")
            def crai = file("${cram}.crai")
            tuple(s, cram, crai)
        }
        .filter { s, cram, crai -> cram.exists() && crai.exists() }
        .set { ch_existing_cram }

    // TODO:: to pass validation the CRAM readgroups must contain all FASTQ readgroups for that sample
    ch_existing_cram
        .set{ ch_validated_cram }

    /* 
        Find and validate any pre-existing GVCFs
    */
    
    ch_sample_names
        .map { s ->
            def gvcf = file("output/results/vcf/gvcf/${s}.g.vcf.gz")
            def tbi = file("${gvcf}.tbi")
            tuple(s, gvcf, tbi)
        }
        .filter { s, gvcf, tbi -> gvcf.exists() && tbi.exists() }
        .set { ch_existing_gvcf }

    // Validate GVCFs
    // TODO:: to pass validation the the comment line must contain all FASTQ readgroups for that sample
    ch_existing_gvcf
        .set{ ch_validated_gvcf }
    
    emit: 
    validated_fastq = ch_validated_fastq
    validated_cram = ch_validated_cram
    validated_gvcf = ch_validated_gvcf

}

