nextflow.enable.types = true

include {
    SamplesheetRow
    LibraryReads
    SampleReads
} from '../types'


def parseSamplesheetRow(row: SamplesheetRow) {
    def sample = row.sample.trim()
    def lib = row.lib?.trim() ?: sample
    def pop = row.pop?.trim() ?: 'unknown'
    def fwd = row.fwd.trim()
    def rev = row.rev?.trim() ?: ''

    /*
     * Determine input type.
     */
    def fwd_is_url =
        fwd ==~ /(?i)^(https?|ftp):\/\/.+/

    def rev_is_url =
        rev &&
        rev ==~ /(?i)^(https?|ftp):\/\/.+/

    def fwd_is_accession =
        fwd ==~ /(?i)^(SRR|ERR|DRR)\d+$/

    def rev_is_accession =
        rev &&
        rev ==~ /(?i)^(SRR|ERR|DRR)\d+$/

    /*
     * Run accession input.
     */
    if (fwd_is_accession) {
        if (rev) {
            error(
                "Run accession '${fwd}' for sample '${sample}' " +
                "must not have a value in the rev column."
            )
        }

        return record(
            sample: sample,
            lib: lib,
            pop: pop,
            source: 'accession',
            input1: fwd,
            input2: '',
            local_reads: []
        )
    }

    /*
     * Remote URL input.
     */
    if (fwd_is_url) {
        if (!rev_is_url) {
            error(
                "URL input for sample '${sample}' requires URLs " +
                "in both fwd and rev columns; found rev='${rev}'."
            )
        }

        return record(
            sample: sample,
            lib: lib,
            pop: pop,
            source: 'url',
            input1: fwd,
            input2: rev,
            local_reads: []
        )
    }

    /*
     * Local FASTQ input.
     */
    if (!rev || rev_is_url || rev_is_accession) {
        error(
            "Local input for sample '${sample}' requires local " +
            "files in both fwd and rev columns; found rev='${rev}'."
        )
    }

    def r1 = file(
        fwd,
        checkIfExists: true
    )

    def r2 = file(
        rev,
        checkIfExists: true
    )

    return record(
        sample: sample,
        lib: lib,
        pop: pop,
        source: 'local',
        input1: r1.name,
        input2: r2.name,
        local_reads: [r1, r2]
    )
}


def groupSampleReads(
    sample: String,
    reads: Bag<LibraryReads>
) {

    def sources = reads
        .collect { read ->
            read.source
        }
        .toSet()

    if (sources.size() != 1) {
        error(
            "Sample '${sample}' contains multiple input source types: " +
            "${sources.join(', ')}."
        )
    }

    def pops = reads
        .collect { read ->
            read.pop
        }
        .toSet()

    if (pops.size() != 1) {
        error(
            "Sample '${sample}' has conflicting population assignments: " +
            "${pops.join(', ')}."
        )
    }

    /*
     * Ensure library ordering is deterministic.
     */
    def sorted_reads = reads.toSorted { read ->
        read.lib
    }

    def source = sources.toList()[0]
    def pop = pops.toList()[0]

    return record(
        sample: sample,
        pop: pop,
        libs: sorted_reads.collect { read ->
            read.lib
        },
        source: source,
        input1s: sorted_reads.collect { read ->
            read.input1
        },
        input2s: sorted_reads.collect { read ->
            read.input2
        },
        local_r1s: source == 'local'
            ? sorted_reads.collect { read ->
                read.local_reads[0]
            }
            : [],
        local_r2s: source == 'local'
            ? sorted_reads.collect { read ->
                read.local_reads[1]
            }
            : []
    )
}
