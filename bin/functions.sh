# This script contains bash functions re-used across processes


#-------------------------------------------------------------------------------
# check_pipeline
#
# Validate the exit status of the most recently executed pipeline.
#
# This function inspects Bash's PIPESTATUS array and returns the first
# non-zero exit code encountered. It is intended for use immediately after
# a pipeline when running with `set -o pipefail`, allowing Nextflow to
# receive the correct exit status and apply retry logic for transient
# failures (e.g. OOM kills).
#
# The function also prints the full PIPESTATUS array to stderr to aid
# debugging of multi-stage pipelines.
#
# Returns:
#   0   All pipeline stages completed successfully.
#   N   Exit code of the first failed stage.
#
# Common exit codes:
#   137  Process killed with SIGKILL (often OOM)
#   143  Process terminated with SIGTERM (walltime/scheduler cancellation)
#   139  Segmentation fault
#
# Notes:
#   - Must be called immediately after the pipeline.
#   - Any intervening command will overwrite PIPESTATUS.
#
# Example:
#   set +e
#   bcftools mpileup ... \
#       | bcftools call ... \
#       | bcftools view ...
#   st=("${PIPESTATUS[@]}")
#   set -e
#   check_pipeline "${st[@]}" || exit $?
#
#-------------------------------------------------------------------------------
check_pipeline() {
    local -a st=("$@")

    echo "PIPESTATUS: ${st[*]}" >&2

    for i in "${!st[@]}"; do
        if (( st[i] != 0 )); then
            echo "Pipeline stage $((i + 1)) failed with exit code ${st[i]}" >&2
            return "${st[i]}"
        fi
    done

    return 0
}

# Download fastq using aria2c
download_fastq_stream_aria2c() {
    local url="$1"
    local threads="${2:-1}"

    aria2c \
        --max-tries=0 \
        --retry-wait=30 \
        --timeout=60 \
        --connect-timeout=60 \
        --file-allocation=none \
        -x "${threads}" \
        -s "${threads}" \
        -o - \
        "${url}"
}

# Download fastq using curl
download_fastq_stream_curl() {
    local url="$1"

    curl \
        --fail \
        --silent \
        --show-error \
        --location \
        --retry 10 \
        --retry-delay 30 \
        --retry-all-errors \
        --connect-timeout 60 \
        --max-time 21600 \
        "${url}"
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
        echo \
            "ERROR: URL did not return gzip data: ${url}" \
            >&2
        echo \
            "ERROR: expected gzip signature 1f8b, received '${magic:-no data}'" \
            >&2
        return 1
    fi

    return 0
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