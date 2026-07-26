process VALIDATE_FASTQ {

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), stdout, emit: status
    
    script:
    """
    #!/usr/bin/env bash
    
    # Returns 0 (true) if gzip is valid AND does NOT have trailing garbage
    check_trailing() {
        local file="\$1"
        local out

        if out=\$(gzip -t -- "\${file}" 2>&1); then
            ! grep -qi 'trailing garbage' <<< "\${out}"
        else
            return 1
        fi
    }

    OK_F=0; OK_R=0; OK_PAIRED=0
    check_trailing "${fastq1}" && OK_F=1
    check_trailing "${fastq2}" && OK_R=1

    # Compare read IDs only after both compressed files pass validation.
    # diff -q exits at the first mismatch.
    if (( OK_F && OK_R )); then
        if diff -q \
            <(seqkit seq -n -i "${fastq1}" | sed -E 's|(/[12])\$||' ) \
            <(seqkit seq -n -i "${fastq2}" | sed -E 's|(/[12])\$||' ) \
            >/dev/null
        then
            OK_PAIRED=1
        fi
    fi

    if (( OK_F && OK_R && OK_PAIRED )); then
        STATUS="PASS"
    else
        STATUS="FAIL"
    fi

    # Extract flow-cell and lane information when R1 is readable.
    if (( OK_F )); then
        READ_HEADER=\$(
            seqkit head -n 1 "${fastq1}" |
                sed -n '1{s|/1\$||;p;}'
        )

        # SRA headers do not contain Illumina flow-cell and lane fields.
        if [[ "\${READ_HEADER}" == @SRR* ]]; then
            FCID="SRA"
            LANE="SRA"
        else
            IFS=':' read -r _ _ FCID LANE _ <<< "\${READ_HEADER}"

            FCID="\${FCID:-UNKNOWN}"
            LANE="\${LANE:-UNKNOWN}"
        fi
    else
        FCID="UNKNOWN"
        LANE="UNKNOWN"
    fi

    # TODO: Detect platform from metadata rather than assuming Illumina.
    PLATFORM="ILLUMINA"

    # Print information to STDOUT to be split into tuples using nextflow logic
    printf '%s\\t%s\\t%s\\t%s\\n' \
        "\${FCID}" \
        "\${LANE}" \
        "\${PLATFORM}" \
        "\${STATUS}"

    """
}