nextflow.enable.types = true
record SamplesheetRow {
    sample: String
    lib: String?
    pop: String?
    fwd: String
    rev: String?
}

record LibraryReads {
    sample: String
    lib: String
    pop: String
    source: String
    input1: String
    input2: String
    local_reads: List<Path>
}

record SampleReads {
    sample: String
    pop: String
    libs: List<String>
    source: String
    input1s: List<String>
    input2s: List<String>
    local_r1s: List<Path>
    local_r2s: List<Path>
}

record SampleCram {
    sample: String
    cram: Path
    crai: Path
}

record CramStats {
    sample: String
    stats: Path
    wgs_coverage: Path
    wgs_per_base: Path
    mean_coverage: String
}

record PerBaseDepth {
    sample: String
    bigwig: Path
}

record PreparedGenome {
    ref_genome: Path
    ref_genome_uri: String
    genome_index_files: Set<Path>
    genome_bed: Path
    reference_masks: Path
    reference_bases: Path
    long_bed: Path
    short_bed: Path
}

record PreparedMito {
    mito_genome: Path
    mito_index_files: Set<Path>
    mito_bed: Path
    shifted_mito_genome: Path
    shifted_mito_index_files: Set<Path>
    shifted_mito_bed: Path
    included_bases: Path
}

record MaskSummary {
    mask_summary_bed: Path
    mask_summary_bed_tbi: Path
    mask_pass_bed: Path
    mask_pass_bed_tbi: Path
    mask_summary: Path
}

record MitoConsensus {
    consensus_fasta: Path
    calls: Path
    qc: Path
}

record PlinkDataset {
    outname: String
    bed: Path
    bim: Path
    fam: Path
}

record VariantChunk {
    interval_id: String
    interval_bed: Path
    interval_index: Path
    variants: Path
    variants_index: Path
}

record VariantFile {
    outname: String
    vcf: Path
    index: Path
}

// TODO - If i move this to VCF can use the existing VariantChunk record instead
record MitoPileup {
    reference: String
    interval_id: String
    samples_tsv: Path
    counts_tsv: Path
}

record SampleGvcf {
    sample: String
    vcf: Path
    index: Path
}

record UnmappedReads {
    sample: String
    r1: Path
    r2: Path
}
