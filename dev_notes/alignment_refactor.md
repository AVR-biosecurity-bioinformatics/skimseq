# Skimseq alignment rework

Implementation plan:
- Remove fastq splitting & validation code, keep bam merging
- Require local fastqs to start
- Implement this with single fastq input
- Update to work with mutliple input fastqs and merge
- Need to update BAM validation

How to deal with read groups
- Still pull from reads?
- Obtain from ENA metadata?


# Handling read groups

Use auxilliary tag in fastq i.e.
@read_name RG:Z:FC001.L001.lib1

Or if thats not maintained through all tools in the pipe, can put it in the read name itself
@FC001.L001.lib1|original_read_name
@FC002.L003.lib2|original_read_name
```

# Function to annotate fastqs with RG's
annotate_fastq() {
    local rg=$1

    awk -v rg="${rg}" '
        NR % 4 == 1 {
            sub(/^@/, "")
            print "@" rg "|" $0
            next
        }
        { print }
    '
}

# Inject the headers before the fifo
{
    stream_rg1_r1 | annotate_fastq "FC001.L001.lib1"
    stream_rg2_r1 | annotate_fastq "FC002.L003.lib2"
} > merged_R1.fifo &

{
    stream_rg1_r2 | annotate_fastq "FC001.L001.lib1"
    stream_rg2_r2 | annotate_fastq "FC002.L003.lib2"
} > merged_R2.fifo &


# Function to inject headers after alignment
inject_sam_readgroups() {
    local rg_header_file=$1

    awk '
        # Read the @RG records from the first input file.
        NR == FNR {
            if ($0 != "") {
                rg_headers[++n_rg] = $0
            }
            next
        }

        # Pass existing SAM header lines through unchanged.
        /^@/ {
            print
            next
        }

        # Immediately before the first alignment, emit all @RG records.
        !headers_injected {
            for (i = 1; i <= n_rg; i++) {
                print rg_headers[i]
            }

            headers_injected = 1
        }

        # Pass alignment records through unchanged.
        {
            print
        }

        END {
            # Handle an unusual header-only SAM stream.
            if (!headers_injected) {
                for (i = 1; i <= n_rg; i++) {
                    print rg_headers[i]
                }
            }
        }
    ' "${rg_header_file}" -
}

| minibwa map \
    -x ${params.minibwa_preset} \
    -k ${params.minibwa_min_seed_length} \
    -c ${params.minibwa_max_seed_occurrence} \
    -t "${aln_threads}" \
    "${ref_genome}" \
    - \
| inject_sam_readgroups readgroups.sam \
| dupblaster \
```



# Then in validate_cram

Could i make this run only if there is not an appropriate cram in work directory?
Or is it better to have a cram input option?

```
def getAlignParamHash(ref_genome) {
    [
        "preset=${params.minibwa_preset}",
        "min_seed=${params.minibwa_min_seed_length}",
        "max_occ=${params.minibwa_max_seed_occurrence}",
        "reference=${ref_genome.baseName}"
    ].join(';').md5()
}
ACTUAL_HASH=$(
    samtools view -H "${cram}" \
    | sed -n 's/^@CO.*PARAM_HASH://p'
)

if [[ "${ACTUAL_HASH}" != "${EXPECTED_HASH}" ]]; then
    STATUS="FAIL"
fi

# COuld check program groups
validate_cram_parameters() {
    local cram=$1
    shift

    local header
    header=$(samtools view -H "${cram}") || return 1

    local expected
    local failed=0

    for expected in "$@"; do
        if grep -Fq -- "${expected}" <<< "${header}"; then
            printf 'PASS: %s\n' "${expected}"
        else
            printf 'FAIL: missing from header: %s\n' "${expected}" >&2
            failed=1
        fi
    done

    return "${failed}"
}
validate_cram_parameters \
    sample.cram \
    'minibwa map' \
    '-x sr' \
    '-k 19' \
    '-c 500' \
    'dupblaster' \
    'samtools sort'

```

# Required conda packages
- conda-forge::aria2=1.37.0
