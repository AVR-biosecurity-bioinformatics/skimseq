# Skimseq alignment rework

Implementation plan:
- Remove fastq splitting & validation code, keep bam merging
- Implement this with single fastq input
- Update to work with mutliple input fastqs and merge
- Need to update BAM validation

How to deal with read groups
- Still pull from reads?
- Obtain from ENA metadata?


# FIFO donwload functions (in a separate file to be sourced)

Could go in bin/functions.sh

```bash
download_fastq_stream() {
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

# CURL VERSION FOR TESTING
download_fastq_stream() {
    local url="$1"

    curl \
        --fail \
        --silent \
        --show-error \
        --location \
        --retry 999999 \
        --retry-delay 30 \
        "${url}"
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
```

# Nextflow process
Download/alignment process

TODO: How to handle multiple input paths, either multiple fastqs from same lane, or multiple remotes with same name
TODO: How to handle bam / read validation in this new format?
```bash
# Source bash functions
source ./${functions}

# Default read sources
FASTQ1="${lib}_R1.fastq.gz"
FASTQ2="${lib}_R2.fastq.gz"

# Handle different fastq sources
if [[ "${source}" == "local" ]]; then
    # Local files
    FASTQ1="${input1}"
    FASTQ2="${input2}"
else
    # Files to be downloaded from online repositories
    case "${source}" in
        url)
            URL1="${input1}"
            URL2="${input2}"
            ;;
        ena)
            read -r URL1 MD5_1 URL2 MD5_2 < <(
                resolve_fastqs "${input1}"
            )
            ;;
        sra)
            read -r URL1 MD5_1 URL2 MD5_2 < <(
                resolve_fastqs "${input1}"
            )
            ;;            
        *)
            echo "ERROR: unsupported input source: ${source}" >&2
            exit 1
            ;;
    esac

    # Make named pipe (consumed by alignment)
    mkfifo "${FASTQ1}" "${FASTQ2}"

    # Download fastq1 as stream
    download_fastq_stream "${URL1}" "${THREADS}" > "${FASTQ1}" &
    PID1=$!
    
    # Download fastq2 as stream
    download_fastq_stream "${URL2}" "${THREADS}" > "${FASTQ2}" &
    PID2=$!

fi

# Create read groups for sample
RG="@RG\\\\tID:\${SAMPLE}.\${LIB}.\${input_id}\\\\tSM:\${SAMPLE}\\\\tLB:\${LIB}\\\\tPU:\${input_id}\\\\tPL:\${platform}"

# Trim adapters, align, mark duplicates, output sorted cram
# In case of corrupted fastq, seqkit sana fixes but pairs may become out of sync
# FASTP should handle this
        
fastp \
    --in1 <(seqkit sana --threads "${THREADS}" "${FASTQ1}") \
    --in2 <(seqkit sana --threads "${THREADS}" "${FASTQ2}") \
    --disable_trim_poly_g \
    --disable_quality_filtering \
    --disable_length_filtering \
    --dont_eval_duplication \
    --stdout \
    --thread ${task.cpus} \
    --json "\${LIB}.\${input_id}.fastp.json" \
    --html "\${LIB}.\${input_id}.fastp.html" \
| minibwa map \
    -x ${params.minibwa_preset} \
    -k ${params.minibwa_min_seed_length} \
    -c ${params.minibwa_max_seed_occurrence} \
    -t "\${ALN_THREADS}" \
    -R '${read_group}' \
    "${ref_genome}" \
    - \
| dupblaster \
| samtools sort \
    -@ ${task.cpus} \
    -O CRAM \
    --reference "${REF}" \
    -o "${lib}.cram"

samtools index -@ ${task.cpus} "${lib}.cram"

# Ensure any background downloads completed successfully
case "${source}" in
    url|ena|sra)
        wait "${PID1}"
        wait "${PID2}"
        ;;
esac

# Clean up FIFOs
rm -f "${FASTQ1}" "${FASTQ2}"

```

# Validation
Can keep md5 or fastq names in the header, i.e

Maybe filtering parameters too? 
```
# Note - add additional filtering parameters etc
def param_string = [
    "preset=${params.minibwa_preset}",
    "min_seed=${params.minibwa_min_seed_length}",
    "max_occ=${params.minibwa_max_seed_occurrence}"
    reference_md5=${ref_md5}
].join(';')

def param_hash = param_string.md5()

def comments = [
    "@CO\tSAMPLE:${sample}",
    "@CO\tLIBRARIES:${rg_list.collect { it[1] }.unique().sort().join(',')}",
    "@CO\tFASTQS:${(fastq1 + fastq2).collect { it.name }.sort().join(',')}",
    "@CO\tPARAM_HASH:${param_hash}",
    "@CO\tPARAMS:${param_string}"

].join('\n')

samtools view -H input.cram > header.sam

cat >> header.sam << 'EOF'
${comments}
EOF

samtools reheader header.sam input.cram > output.cram
samtools index output.cram

```

# Then in validate_cram
````
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
```

# TESTING

Below works!
```bash
ml SeqKit/2.8.2
ml fastp/0.23.4-GCC-13.3.0
ml bwa-mem2/2.2.1-GCC-13.3.0
ml SAMtools/1.23.1-GCC-13.3.0

accession=SRR21855497

rm $FASTQ1 $FASTQ2
FASTQ1="${accession}_R1.fastq.gz"
FASTQ2="${accession}_R2.fastq.gz"

ref_genome=/group/referencedata/mspd-db/genomes/insect/halyomorpha_halys/GCF_000696795.3_Hhal_1.1_genomic.fna
THREADS=2
 read -r URL1 MD5_1 URL2 MD5_2 < <(
     resolve_fastqs "${accession}"
)
# Make named pipe (consumed by alignment)
mkfifo "${FASTQ1}" "${FASTQ2}"

# Download fastq1 as stream
download_fastq_stream "${URL1}" 1 > "${FASTQ1}" &
PID1=$!
    
# Download fastq2 as stream
download_fastq_stream "${URL2}" 1 > "${FASTQ2}" &
PID2=$!


fastp \
    --in1 <(seqkit sana --threads 1 "${FASTQ1}") \
    --in2 <(seqkit sana --threads 1 "${FASTQ2}") \
    --disable_trim_poly_g \
    --disable_quality_filtering \
    --disable_length_filtering \
    --dont_eval_duplication \
    --stdout \
    --thread 1 \
    --json "\${LIB}.\${input_id}.fastp.json" \
    --html "\${LIB}.\${input_id}.fastp.html" \
| bwa-mem2 mem \
    -t 4 \
    "${ref_genome}" \
    - \
| samtools view \
    -@ 1 \
    -O CRAM \
    --reference "${ref_genome}" \
    -o "${accession}.cram"

wait "${PID1}"
wait "${PID2}"
```


# Required conda packages
- conda-forge::aria2=1.37.0



# OLD FUNCTIONS PRE FIFO
functions.sh
```bash

# Downloads one file using aria2c with automatic retries and resumption
# Usage: download_fastq URL OUTPUT [MD5] [THREADS]
# aria2c retries/resumes until transfer succeeds So dont need to retry failure on download
# Dependencies: aria2c
download_fastq() {
    local url="$1"
    local output="$2"
    local md5="${3:-}"
    local threads="${4:-1}"
    local checksum_args=()

    (( threads > 0 )) || {
        echo "ERROR: threads must be greater than zero" >&2
        return 1
    }

    if [[ -n "${md5}" ]]; then
        checksum_args+=(
            --check-integrity=true
            "--checksum=md5=${md5}"
        )
    fi

    aria2c \
        -c \
        -x "${threads}" \
        -s "${threads}" \
        --retry-wait=30 \
        --max-tries=10 \
        --timeout=60 \
        --connect-timeout=60 \
        --file-allocation=none \
        "${checksum_args[@]}" \
        -o "${output}" \
        "${url}"

    [[ -s "${output}" ]] || {
        echo "ERROR: downloaded file is missing or empty: ${output}" >&2
        return 1
    }
}
# Downloads an NCBI SRA run using prefetch, validates, then converts to gzipped fastq
# Usage: download_sra_fastqs ACCESSION [OUTPUT_PREFIX] [THREADS] [TMPDIR]
# Dependencies: prefetch, vdb-validate, fasterq-dump, pigz
download_sra_fastqs() {
    local accession="$1"
    local output_prefix="${2:-${accession}}"
    local threads="${3:-1}"
    local tmpdir="${4:-${TMPDIR:-.}/fasterq.${accession}}"

    # Check inputs 
    [[ "${accession}" =~ ^SRR[0-9]+$ ]] || {
        echo "ERROR: invalid SRA run accession: ${accession}" >&2
        return 1
    }

    [[ "${threads}" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: threads must be a positive integer: ${threads}" >&2
        return 1
    }

    # Download SRA format file
    mkdir -p "${tmpdir}"
    prefetch --max-size u --output-directory "${accession}" "${accession}"

    # Validate downloaded file
    vdb-validate "${accession}"

    # Convert to fastq
    fasterq-dump \
        --split-files \
        --threads "${threads}" \
        --temp "${tmpdir}" \
        "${accession}"

    # Gzip fastqs
    pigz --processes "${threads}" "${accession}_1.fastq" "${accession}_2.fastq"

    # Check fastq.gz are present
    [[ -s "${accession}_1.fastq.gz" ]] || { echo "ERROR: compressed R1 FASTQ is missing or empty" >&2 return 1 }
    [[ -s "${accession}_2.fastq.gz" ]] || { echo "ERROR: compressed R2 FASTQ is missing or empty" >&2 return 1 }

    # Rename fastqs
    if [[ "${output_prefix}" != "${accession}" ]]; then
        mv -f "${accession}_1.fastq.gz" "${output_prefix}_1.fastq.gz"
        mv -f "${accession}_2.fastq.gz" "${output_prefix}_2.fastq.gz"
    fi

    # Clean up temporary directory
    rm -rf "${tmpdir}"
}

# Resolves an ENA or DDBJ accession through ENA Portal API, extracts paired FASTQ URLs  MD5 checksums, then downloads both using download_fastq().
# Usage: download_ena_fastqs ACCESSION [OUTPUT_PREFIX] [THREADS]
# Dependencies: curl, awk, download_fastq function
download_ena_fastqs() {
    local accession="$1"
    local output_prefix="${2:-${accession}}"
    local threads="${3:-1}"

    local api_url
    local metadata
    local fastq_urls
    local fastq_md5s

    local -a urls=()
    local -a md5s=()

    [[ "${accession}" =~ ^(ERR|DRR)[0-9]+$ ]] || {
        echo "ERROR: invalid ENA/DDBJ run accession: ${accession}" >&2
        return 1
    }

    [[ "${threads}" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: threads must be a positive integer: ${threads}" >&2
        return 1
    }

    # Retrieve metadata from api
    api_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    metadata=$(
        curl \
            --fail \
            --silent \
            --show-error \
            --location \
            --retry 5 \
            --retry-delay 10 \
            "${api_url}"
    )

    [[ -n "${metadata}" ]] || {
        echo "ERROR: ENA returned no metadata for ${accession}" >&2
        return 1
    }

    IFS=$'\t' read -r fastq_urls fastq_md5s < <(
        printf '%s\n' "${metadata}" |
            awk -F '\t' 'NR == 2 { print $1 "\t" $2 }'
    )

    [[ -n "${fastq_urls}" ]] || {
        echo "ERROR: ENA returned no FASTQ URLs for ${accession}" >&2
        return 1
    }

    IFS=';' read -r -a urls <<< "${fastq_urls}"

    if [[ -n "${fastq_md5s}" ]]; then
        IFS=';' read -r -a md5s <<< "${fastq_md5s}"
    fi

    if (( ${#urls[@]} != 2 )); then
        echo "ERROR: expected exactly two FASTQ files for ${accession}, found ${#urls[@]}" >&2
        return 1
    fi

    if [[ -n "${fastq_md5s}" ]] && (( ${#md5s[@]} != ${#urls[@]} )); then
        echo "ERROR: FASTQ URL and MD5 counts differ for ${accession}" >&2
        return 1
    fi

    [[ "${urls[0]}" == *://* ]] || urls[0]="ftp://${urls[0]}"
    [[ "${urls[1]}" == *://* ]] || urls[1]="ftp://${urls[1]}"

    download_fastq \
        "${urls[0]}" \
        "${output_prefix}_1.fastq.gz" \
        "${md5s-}" \
        "${threads}"

    download_fastq \
        "${urls[1]}" \
        "${output_prefix}_2.fastq.gz" \
        "${md5s-}" \
        "${threads}"

    [[ -s "${output_prefix}_1.fastq.gz" ]] || {
        echo "ERROR: ENA R1 download is missing or empty" >&2
        return 1
    }

    [[ -s "${output_prefix}_2.fastq.gz" ]] || {
        echo "ERROR: ENA R2 download is missing or empty" >&2
        return 1
    }
}

# Validation helpers

# Extracts FASTQ record IDs and removes trailing /1 or /2 mate suffixes 
# Usage: normalise_ids FASTQ [THREADS]
# Dependencies: seqkit, sed
normalise_ids() {
    local fastq="$1"
    local threads="${2:-1}"

    seqkit seq \
        --threads "${threads}" \
        -n \
        -i \
        "${fastq}" |
        sed -E 's|(/[12])$||'
}

# Checks that two FASTQ.gz files are non-empty and contain valid gzip streams
# Usage: validate_structure R1 R2 [THREADS]
# Dependencies: gzip, seqkit
validate_gzip() {
    local fq1="$1"
    local fq2="$2"

    [[ -s "${fq1}" && -s "${fq2}" ]] || return 1

    gzip -t -- "${fq1}" || return 1
    gzip -t -- "${fq2}" || return 1
}

# Verifies that R1 and R2 contain the same normalised read IDs in the same order. Stops at first mismatch
# Usage: validate_pairing R1 R2 [THREADS]
# Dependencies: cmp, normalise_ids, seqkit, sed
validate_pairing() {
    local fq1="$1"
    local fq2="$2"
    local threads="${3:-1}"

    cmp -s \
        <(normalise_ids "${fq1}" "${threads}") \
        <(normalise_ids "${fq2}" "${threads}")
}

# Performs complete paired-FASTQ validation by combining structural and read-pair synchronisation checks.
# Usage: validate_fastqs R1 R2 [THREADS]
# Dependencies: validate_structure, validate_pairing
validate_fastqs() {
    local fq1="$1"
    local fq2="$2"
    local threads="${3:-1}"

    validate_gzip "${fq1}" "${fq2}" || return 1
    validate_pairing "${fq1}" "${fq2}" "${threads}" || return 1
}


# Repair malformed FASTQ records and restore pairing
# Usage: repair_fastqs R1 R2 [THREADS] [TMPDIR]
repair_fastqs() {

    local fq1="$1"
    local fq2="$2"
    local threads="${3:-1}"
    local tmpdir="repaired"
    local rescued1="rescued_1.fq.gz"
    local rescued2="rescued_2.fq.gz"
    local repaired1="${tmpdir}/${rescued1}"
    local repaired2="${tmpdir}/${rescued2}"

    # Sanitise fastq records
    seqkit sana --threads "${threads}" ${fq1} -o ${rescued1}
    seqkit sana --threads "${threads}" ${fq2} -o ${rescued2}
    
    [[ -s ${rescued1} ]] || { echo "ERROR: no R1 reads recovered from ${fq1}" >&2 return 1 }
    [[ -s ${rescued2} ]] || { echo "ERROR: no R2 reads recovered from ${fq2}" >&2 return 1 }
    
    # re-pair sanitised fastqs
    seqkit pair --threads $"${threads}" -1 ${rescued1} -2 ${rescued2} -O ${tmpdir}
    [[ -s "${repaired1}" ]] || { echo "ERROR: pair repair produced no R1 output" >&2 return 1 }
    [[ -s "${repaired2}" ]] || { echo "ERROR: pair repair produced no R2 output" >&2 return 1 }

    # Rename outputs and clean up
    mv -f "${repaired1}" "${fq1}"
    mv -f "${repaired2}" "${fq2}"
    rm -f "${rescued1}" "${rescued2}"
    rmdir "${tmpdir}"
 }
```