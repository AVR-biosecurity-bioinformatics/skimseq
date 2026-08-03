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

    main:
    startupMessage()

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

    publish:

    mask_summary    = SKIMSEQ.out.mask_summary
    mask_summary_bed = SKIMSEQ.out.mask_summary_bed
    mask_pass_bed = SKIMSEQ.out.mask_pass_bed

    cram            = SKIMSEQ.out.cram
    perbase         = SKIMSEQ.out.perbase
    mito_fasta      = SKIMSEQ.out.mito_fasta

    unfiltered_vcf = SKIMSEQ.out.unfiltered_vcf
    gvcf = SKIMSEQ.out.gvcf
    final_vcf = SKIMSEQ.out.final_vcf

    beagle_gl       = SKIMSEQ.out.beagle_gl
    plink           = SKIMSEQ.out.plink
    pca             = SKIMSEQ.out.pca
    relationship    = SKIMSEQ.out.relationship
    king            = SKIMSEQ.out.king
    distance        = SKIMSEQ.out.distance
    ordination_plot = SKIMSEQ.out.ordination_plot
    pca_plot        = SKIMSEQ.out.pca_plot
    tree_plot       = SKIMSEQ.out.tree_plot
    popmap          = SKIMSEQ.out.popmap

    // QC outputs
    sample_filter_plots = SKIMSEQ.out.sample_filter_plots
    sample_missing_tsv = SKIMSEQ.out.sample_missing_tsv
    site_filter_plots = SKIMSEQ.out.site_filter_plots
    cram_stats      = SKIMSEQ.out.cram_stats
    cram_plots      = SKIMSEQ.out.cram_plots
    vcf_stats       = SKIMSEQ.out.vcf_stats
    multiqc_report  = SKIMSEQ.out.multiqc_report
    multiqc_plots   = SKIMSEQ.out.multiqc_plots
    multiqc_data    = SKIMSEQ.out.multiqc_data

}

output {
    mask_summary {
        path 'qc'
    }
    mask_summary_bed {
        path 'qc'
    }
    mask_pass_bed {
        path 'qc'
    }
    cram {
        enabled params.output_cram
        path params.cram_store.startsWith("${workflow.outputDir}/")
            ? params.cram_store - "${workflow.outputDir}/"
            : params.cram_store
    }
    gvcf {
        enabled params.output_gvcf
        path params.gvcf_store.startsWith("${workflow.outputDir}/")
            ? params.gvcf_store - "${workflow.outputDir}/"
            : params.gvcf_store
    }
    mito_fasta {
        path 'mito'
    }
    unfiltered_vcf {
        path 'vcf/unfiltered'
    }
    final_vcf {
        path 'vcf/filtered'
    }

    beagle_gl {
        enabled params.output_beagle_gl
        path 'beagle'
    }
    plink {
        path 'plink'
    }
    pca {
        path 'plink'
    }
    relationship {
        path 'plink'
    }
    king {
        path 'plink'
    }
    distance {
        path 'distmat'
    }
    ordination_plot {
        path 'visualisation/ordination'
    }
    pca_plot {
        path 'visualisation/ordination'
    }
    tree_plot {
        path 'visualisation/trees'
    }
    popmap {
        path 'metadata'
    }
    perbase {
        path 'qc/cram_stats/perbase_depth'
    }
    cram_stats {
        path 'qc/cram_stats/cram_qc_data'
    }
    cram_plots {
        path 'qc/cram_stats/cram_qc_plots'
    }
    vcf_stats {
        path 'qc/vcf_stats'
    }
    sample_filter_plots {
        path 'qc'
    }   
    sample_missing_tsv {
        path 'qc'
    }  
    site_filter_plots {
        path 'qc'
    }  
    multiqc_report {
        path 'qc'
    }   
    multiqc_plots {
        path 'qc'
    }    
    multiqc_data {
        path 'qc'
    }    
}