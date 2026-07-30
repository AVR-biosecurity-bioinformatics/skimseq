#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AVR-biosecurity-bioinformatics/skimseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/AVR-biosecurity-bioinformatics/skimseq
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include functions from nf-schema

include {
    paramsHelp
    validateParameters
    paramsSummaryLog
    samplesheetToList
} from 'plugin/nf-schema'

///// define functions 

def startupMessage() {
    log.info pipelineHeader()
    log.info "~~~ skimseq: A bioinformatics pipeline for processing low or variable coverage whole genome sequencing data ~~~"
    log.info " "
}

def pipelineHeader() {
    return """                                                                          
                       __    __                            
                .-----|  |--|__.--------.-----.-----.-----.
                |__ --|    <|  |        |__ --|  -__|  _  |
                |_____|__|__|__|__|__|__|_____|_____|__   |
                                                       |__|
    """.stripIndent()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// workflows
include { SKIMSEQ                                                   } from './nextflow/workflows/skimseq'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

///// run implicit workflow
workflow {

    startupMessage()
    

    // Print summary of supplied parameters (that differ from defaults)
    log.info paramsSummaryLog(workflow)

    // validate inpiut params
    validateParameters()
    log.info paramsSummaryLog(workflow)

    workflow.onComplete = {
        if ( workflow.success ) {
        log.info "[$workflow.complete] >> Pipeline finished SUCCESSFULLY after $workflow.duration"
        } else {
        log.info "[$workflow.complete] >> Pipeline finished with ERRORS after $workflow.duration"
        }
    }

    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
        //    log.info startupMessage()
        log.info paramsHelp("nextflow run AVR-biosecurity-bioinformatics/skimseq") // TODO: add typical commands for pipeline
        exit 0
    }


    ///// run main skimseq workflow
    SKIMSEQ ()

}