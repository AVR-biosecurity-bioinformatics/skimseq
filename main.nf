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

    // Hack to enforce cram_store and gvcf_store to be within results directory
    // TODO: Replace directory-based reuse with CRAM/GVCF manifests.
    // Workflow outputs currently require publication beneath workflow.outputDir.
    def outdir = new File(workflow.outputDir.toString()).canonicalPath

    if( params.output_cram ) {
        def cramStore = new File(params.cram_store).canonicalPath

        if( !cramStore.startsWith(outdir) ) {
            error """
            params.cram_store must be located within workflow.outputDir

            workflow.outputDir = ${outdir}
            cram_store         = ${cramStore}

            External CRAM output locations are currently unsupported.
            """
        }
    }

    if( params.output_gvcf ) {
        def gvcfStore = new File(params.gvcf_store).canonicalPath

        if( !gvcfStore.startsWith(outdir) ) {
            error """
            params.gvcf_store must be located within workflow.outputDir

            workflow.outputDir = ${outdir}
            gvcf_store         = ${gvcfStore}

            External GVCF output locations are currently unsupported.
            """
        }
    }

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
    newick_tree     = SKIMSEQ.out.newick_tree
    popmap          = SKIMSEQ.out.popmap

    // QC outputs
    sample_filter_plots = SKIMSEQ.out.sample_filter_plots
    sample_missing_tsv = SKIMSEQ.out.sample_missing_tsv
    site_filter_plots = SKIMSEQ.out.site_filter_plots
    cram_stats      = SKIMSEQ.out.cram_stats
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
    // NOTE: For now cram store and gvcf store have to be within results directory
    cram {
        enabled params.output_cram
        path new File(params.cram_store).name
    }
    gvcf {
        enabled params.output_gvcf
        path new File(params.gvcf_store).name
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
    newick_tree {
        path 'visualisation/trees'
    }
    popmap {
        path 'metadata'
    }
    perbase {
        enabled params.output_perbase_depth
        path 'qc/cram_stats/perbase_depth'
    }
    cram_stats {
        path 'qc/cram_stats/cram_qc_data'
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