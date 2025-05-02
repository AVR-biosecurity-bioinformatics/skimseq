process FASTP {
    def process_name = "fastp"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '4.GB'
    cpus 4
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "fastp/0.23.4-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    tuple val(rf_quality), val(rf_length), val(rf_n_bases), val(rf_trim_polyg), val(rf_cut_right), val(rf_cut_window_size), val(rf_cut_mean_quality), val(rf_lc_filter), val(rf_lc_threshold), val(rf_correction), val(rf_overlap_length), val(rf_overlap_diff), val(rf_overlap_diff_pc), val(rf_custom_flags)

    output: 
    tuple val(sample), path("*.trimmed.R1.fastq"), path("*.trimmed.R2.fastq"), path("*.json"),     emit: fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${fastq1} \
        ${fastq2} \
        ${rf_quality} \
        ${rf_length} \
        ${rf_n_bases} \
        ${rf_trim_polyg} \
        ${rf_cut_right} \
        ${rf_cut_window_size} \
        ${rf_cut_mean_quality} \
        ${rf_lc_filter} \
        ${rf_lc_threshold} \
        ${rf_correction} \
        ${rf_overlap_length} \
        ${rf_overlap_diff} \
        ${rf_overlap_diff_pc} \
        "${rf_custom_flags}"

    """
}