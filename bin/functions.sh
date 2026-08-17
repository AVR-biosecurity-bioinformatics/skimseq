#!/bin/bash

# This script contains bash functions re-used across processes


#-------------------------------------------------------------------------------
# check_pipeline
#
# Validate the exit status of the most recently executed pipeline.
# Prioritises substantive errors before sigpipe
check_pipeline() {
    local -a st=("$@")
    local i

    echo "PIPESTATUS: ${st[*]}" >&2

    # Prioritise substantive failures over SIGPIPE.
    # An upstream 141 is commonly caused by a downstream process failing.
    for i in "${!st[@]}"; do
        if (( st[i] != 0 && st[i] != 141 )); then
            echo \
                "Pipeline stage $((i + 1)) failed with exit code ${st[i]}" \
                >&2
            return "${st[i]}"
        fi
    done

    # Return SIGPIPE only if it is the only type of failure.
    for i in "${!st[@]}"; do
        if (( st[i] == 141 )); then
            echo \
                "Pipeline stage $((i + 1)) failed with SIGPIPE (141)" \
                >&2
            return 141
        fi
    done

    return 0
}

# Download fastq using HydraStream
download_fastq_stream_hs() {
    local url="$1"
    local expected_md5="${2:-}"
    local threads="${3:-4}"
    local log_dir

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "ERROR: HydraStream threads must be a positive integer; " \
            "received '${threads}'" \
            >&2
        return 2
    fi

    url="${url/#ftp:\/\//https:\/\/}"

    log_dir=$(
        mktemp -d \
            "${TMPDIR:-.}/hydrastream.XXXXXX"
    ) || return 1

    local -a args=(
        "${url}"
        --threads "${threads}"
        --stream
        --no-ui
        --output "${log_dir}"
    )

    if [[ -n "${expected_md5}" ]]; then
        args+=(
            --typehash md5
            --checksum "${expected_md5}"
        )
    fi

    hs "${args[@]}"
    local status=$?

    rm -rf "${log_dir}"

    return "${status}"
}

# Streaming helper
stream_fastq() {
    local input="$1"
    local rg_id="$2"
    local expected_md5="${3:-}"
    local download_threads="${4:-4}"

    {
        case "${STREAM_TYPE}" in
            local)
                salvage_fastq_stream "${input}"
                ;;

            remote)
                download_fastq_stream_hs \
                    "${input}" \
                    "${expected_md5}" \
                    "${download_threads}" |
                    gzip -dc |
                    seqkit sana --threads 1 -
                ;;

            *)
                echo "ERROR: unsupported stream type '${STREAM_TYPE}'" >&2
                return 2
                ;;
        esac
    } |
        annotate_fastq "${rg_id}"
}

salvage_fastq_stream() {
    local input="$1"
    local -a status

    set +e

    gzip -dc -- "${input}" |
        seqkit sana --threads 1 -

    status=("${PIPESTATUS[@]}")

    set -e

    # SeqKit itself must complete successfully.
    if (( status[1] != 0 )); then
        echo "ERROR: FASTQ sanitisation failed for '${input}'" >&2
        return "${status[1]}"
    fi

    # A gzip failure means the stream was partially recovered.
    if (( status[0] != 0 )); then
        echo "WARNING: recovered reads from corrupted gzip file '${input}'" >&2
    fi

    return 0
}

# Validate that a file from URL begins with the gzip magic bytes 1f 8b.
# The request is limited to the first two bytes. Do not use
# --retry-all-errors here because closing a short validation stream can
# otherwise cause curl write errors to be retried.
validate_gzip_url() {
    local url="$1"
    local magic

    if [[ -z "${url}" ]]; then
        echo "ERROR: validate_gzip_url received an empty URL" >&2
        return 2
    fi

    echo "Validating remote gzip stream: ${url}" >&2

    magic=$(
        curl \
            --fail \
            --silent \
            --show-error \
            --location \
            --range 0-1 \
            --retry 3 \
            --retry-delay 5 \
            --connect-timeout 30 \
            --max-time 60 \
            "${url}" \
        | od -An -N2 -tx1 \
        | tr -d '[:space:]'
    )

    if [[ "${magic}" != "1f8b" ]]; then
        echo "ERROR: URL did not return gzip data: ${url}" >&2
        echo "ERROR: expected gzip signature 1f8b, received '${magic:-no data}'" >&2
        return 1
    fi

    return 0
}

# Parse a shortread header
parse_shortread_header() {
    local read_header="$1"
    local -a fields

    # Retain only the first whitespace-delimited component.
    read_header="${read_header#@}"
    read_header="${read_header%%[[:space:]]*}"
    read_header="${read_header%/1}"
    read_header="${read_header%/2}"

    IFS=':' read -r -a fields <<< "${read_header}"

    # Expected structure:
    # instrument:run:flowcell:lane:tile:x:y
    if (( ${#fields[@]} < 4 )); then
        echo \
            "ERROR: header has fewer than four colon-delimited fields: " \
            "'${read_header}'" \
            >&2
        return 1
    fi

    local fcid="${fields[2]}"
    local lane="${fields[3]}"

    if [[ -z "${fcid}" || ! "${lane}" =~ ^[0-9]+$ ]]; then
        echo \
            "ERROR: could not parse flowcell and lane from FASTQ header " \
            "'${read_header}'" \
            >&2
        return 1
    fi

    printf '%s %s\n' "${fcid}" "${lane}"
}


# Get local flowcell and lane
get_local_flowcell_lane() {
    local fastq="$1"
    local read_header

    read_header=$(
        gzip -dc -- "${fastq}" |
            head -n 1
    )

    parse_shortread_header "${read_header}"
}

# Get remote flowcell and lane
get_remote_flowcell_lane() {
    local url="$1"
    local read_header
    local qname
    local accession

    read_header=$(
        set +o pipefail

        curl \
            --location \
            --fail \
            --silent \
            --show-error \
            "${url}" 2>/dev/null |
        gzip -dc 2>/dev/null |
        head -n 1
    )

    if [[ -z "${read_header}" ]]; then
        echo "ERROR: could not read FASTQ header from '${url}'" >&2
        return 1
    fi

    qname="${read_header#@}"
    qname="${qname%%[[:space:]]*}"

    # ENA/SRA archive header, e.g. SRR13005336.1
    if [[ "${qname}" =~ ^((SRR|ERR|DRR)[0-9]+)\.[0-9]+$ ]]; then
        accession="${BASH_REMATCH[1]}"
        printf '%s %s\n' "${accession}" "1"
        return 0
    fi

    parse_shortread_header "${read_header}"
}

# Joint flowcell lane parsing function
get_flowcell_lane() {
    local input="$1"
    local stream_type="$2"
    local flowcell_lane
    local fcid
    local lane

    case "${stream_type}" in
        local)
            flowcell_lane=$(
                get_local_flowcell_lane "${input}"
            ) || return 1
            ;;
        remote)
            flowcell_lane=$(
                get_remote_flowcell_lane "${input}"
            ) || return 1
            ;;
        *)
            echo "ERROR: unsupported stream type '${stream_type}' for '${input}'" >&2
            return 1
            ;;
    esac

    read -r fcid lane <<< "${flowcell_lane}"

    if [[ -z "${fcid}" || -z "${lane}" ]]; then
        echo "ERROR: could not determine flowcell/lane for '${input}'" >&2
        return 1
    fi

    printf '%s %s\n' "${fcid}" "${lane}"
}

# Accepts ENA or SRA accession - resolves to ENA
resolve_fastqs() {
    local accession="$1"

    local metadata
    local fastq_urls
    local fastq_md5s

    local -a urls=()
    local -a md5s=()

    [[ "${accession}" =~ ^(SRR|ERR|DRR)[0-9]+$ ]] || {
        echo "ERROR: unsupported accession: ${accession}" >&2
        return 1
    }

    metadata=$(
        curl \
            --fail \
            --silent \
            --show-error \
            --location \
            --retry 5 \
            --retry-delay 10 \
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    )

    read -r accession_returned fastq_urls fastq_md5s <<< "$(
        printf '%s\n' "${metadata}" |
            awk 'NR==2'
    )"

    [[ -n "${fastq_urls}" ]] || {
        echo "ERROR: ENA returned no FASTQ URLs for ${accession}" >&2
        return 1
    }

    IFS=';' read -r -a urls <<< "${fastq_urls}"

    (( ${#urls[@]} == 2 )) || {
        echo "ERROR: expected 2 FASTQ files for ${accession}, found ${#urls[@]}" >&2
        return 1
    }

    if [[ -n "${fastq_md5s}" ]]; then
        IFS=';' read -r -a md5s <<< "${fastq_md5s}"

        (( ${#md5s[@]} == 2 )) || {
            echo "ERROR: expected 2 MD5 values for ${accession}, found ${#md5s[@]}" >&2
            return 1
        }
    fi

    [[ "${urls[0]}" == *://* ]] || urls[0]="ftp://${urls[0]}"
    [[ "${urls[1]}" == *://* ]] || urls[1]="ftp://${urls[1]}"

    printf '%s\t%s\t%s\t%s\n' \
        "${urls[0]}" \
        "${md5s-}" \
        "${urls[1]}" \
        "${md5s-}"
}

# Annotate fastq read names with readgroups
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


# Inject the headers into samtools readgroups
# This can include CO and PG lines
inject_sam_readgroups() {
    local injected_header_file="$1"

    awk '
        BEGIN {
            FS = OFS = "\t"
        }

        # Read records that will be injected into the SAM header.
        NR == FNR {
            if ($0 == "") {
                next
            }

            injected_headers[++n_headers] = $0

            # Only @RG records define valid read-group identifiers.
            if ($1 == "@RG") {
                rg_id = ""

                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^ID:/) {
                        rg_id = substr($i, 4)
                        break
                    }
                }

                if (rg_id == "") {
                    print \
                        "ERROR: @RG record is missing an ID field: " $0 \
                        > "/dev/stderr"
                    exit 1
                }

                valid_rg[rg_id] = 1
                n_rg++
            }

            next
        }

        # Pass the existing SAM header through unchanged.
        /^@/ {
            print
            next
        }

        # Add injected header records immediately before the first alignment.
        !headers_injected {
            if (n_rg == 0) {
                print \
                    "ERROR: auxiliary header contains no @RG records" \
                    > "/dev/stderr"
                exit 1
            }

            for (i = 1; i <= n_headers; i++) {
                print injected_headers[i]
            }

            headers_injected = 1
        }

        {
            # Expected QNAME:
            #
            #   read_group_id|original_read_name
            separator = index($1, "|")

            if (separator == 0) {
                print \
                    "ERROR: alignment QNAME does not contain an RG prefix: " \
                    $1 \
                    > "/dev/stderr"
                exit 1
            }

            rg_id = substr($1, 1, separator - 1)
            original_qname = substr($1, separator + 1)

            if (!(rg_id in valid_rg)) {
                print \
                    "ERROR: QNAME references undefined read group: " \
                    rg_id \
                    > "/dev/stderr"
                exit 1
            }

            if (original_qname == "") {
                print \
                    "ERROR: empty QNAME after removing RG prefix: " \
                    $1 \
                    > "/dev/stderr"
                exit 1
            }

            $1 = original_qname

            # Rebuild the alignment without any pre-existing RG tag.
            output = ""

            for (i = 1; i <= NF; i++) {
                if (i > 11 && $i ~ /^RG:Z:/) {
                    continue
                }

                output = output == "" ? $i : output OFS $i
            }

            print output, "RG:Z:" rg_id
        }

        END {
            # Handle the unusual case of a header-only SAM stream.
            if (!headers_injected) {
                for (i = 1; i <= n_headers; i++) {
                    print injected_headers[i]
                }
            }
        }
    ' "${injected_header_file}" -
}
