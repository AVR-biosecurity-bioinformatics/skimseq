## Phylloxera sample tests
### subsample reads for testing
```
cp /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/*.fastq.gz test

module load SeqKit

# sample 1 (G1)
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R1_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R2_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R2_001.fastq.gz -s1

# sample 2 (G4)
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv18g4vaitc8152/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R1_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv18g4vaitc8152/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R2_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R2_001.fastq.gz -s1



```

### phylloxera genome and mito genome from here: 
`https://bipaa.genouest.org/sp/daktulosphaira_vitifoliae/download/genome/v3.1/`

### test commands

```
export NXF_VER=23.05.0-edge

nextflow run . -profile basc_modules,debug

nextflow run . -profile basc_modules,debug --mito_genome test/Dv_mitochondrial_genome.fa --ref_genome test/Dv_genome_V3.1.fa

# downsampled files 
nextflow run . -profile basc_modules,debug \
    --mito_genome test/Dv_mitochondrial_genome.fa \
    --ref_genome test/Dv_genome_V3.1.fa \
    --samplesheet test/samplesheet.csv

# full files
nextflow run . -profile basc_modules,debug \
    --mito_genome test/Dv_mitochondrial_genome.fa \
    --ref_genome test/Dv_genome_V3.1.fa \
    --samplesheet test/samplesheet_full.csv

```

# Create Qfly test datasets
This test data set uses a small segment of Qfly chromosome 1: CM028320.1:50000-99999, and CM028321.1:50000-59999

```
ml Miniconda3/24.7.1-0
conda create \
    --name rasusa \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority \
    rasusa=5.1.0 \
    minibwa=0.6 \
    samtools=1.24

conda activate rasusa

BAM_DIR="/group/pathogens/IAWS/Personal/Alexp/skimseq_qfly/output/results/cram"
REF_GENOME="/group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/GCA_016617805.2_CSIRO_BtryS06_freeze2_genomic_withmito.fna"
MITO_FASTA="/group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/mitogenome/HQ130030.1_Bactrocera_tryoni_mitochondrion.fa"
OUTPUT_DIR="test_data/qfly"
mkdir -p $OUTPUT_DIR

MITO_CONTIG="HQ130030.1"
SAMPLES=(
    "EM6"
    "EM3"
    "F3"
    "F2xM12-F1"
)

REGIONS=(
    "CM028320.1:50000-99999"
    "CM028321.1:50000-59999"
    "${MITO_CONTIG}"
)

# Build cut down reference
TEST_REFERENCE="${OUTPUT_DIR}/test_qfly_genome.fa"

samtools faidx \
    "$REF_GENOME" \
    "CM028320.1:50000-99999" \
    "CM028321.1:50000-59999" \
    "$MITO_CONTIG" |
    sed 's/:.*$//' \
    > "$TEST_REFERENCE"

samtools faidx "$TEST_REFERENCE"
minibwa index "$TEST_REFERENCE"

for sample in "${SAMPLES[@]}"; do
    cram="${BAM_DIR}/${sample}.cram"

    recruited_r1="${OUTPUT_DIR}/${sample}.recruited_R1.fastq.gz"
    recruited_r2="${OUTPUT_DIR}/${sample}.recruited_R2.fastq.gz"
    recruited_singleton="${OUTPUT_DIR}/${sample}.recruited_singleton.fastq.gz"

    recruited_bam="${OUTPUT_DIR}/${sample}.recruited.bam"
    subsampled_bam="${OUTPUT_DIR}/${sample}.recruited.5x.bam"

    r1="${OUTPUT_DIR}/${sample}_subset_R1.fastq.gz"
    r2="${OUTPUT_DIR}/${sample}_subset_R2.fastq.gz"
    singleton="${OUTPUT_DIR}/${sample}_subset_singleton.fastq.gz"
    category0="${OUTPUT_DIR}/${sample}_subset_cat0.fastq.gz"
    echo "[INFO] Recruiting test reads for $sample" >&2

    ###########################################################################
    # Recruit complete templates from the original CRAM and realign
    ###########################################################################

    samtools view \
        --threads 4 \
        --fetch-pairs \
        --reference "$REF_GENOME" \
        -f 0x1 \
        -F 0xF00 \
        -u \
        "$cram" \
        "${REGIONS[@]}" |
        samtools collate \
            --threads 4 \
            -Ou \
            - |
        samtools fastq \
            --threads 4 \
            -n \
            -0 /dev/null \
            -s /dev/null \
            - |
        minibwa mem \
            -t 4 \
            -p \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
            "$TEST_REFERENCE" \
            - |
        samtools sort \
            --threads 4 \
            -O BAM \
            -o "$recruited_bam" \
            -

    samtools index \
        --threads 4 \
        "$recruited_bam"

    ###########################################################################
    # Subsample every reduced-reference contig to approximately 30x
    ###########################################################################

    rasusa aln \
        --coverage 30 \
        --seed 42 \
        --output "$subsampled_bam" \
        "$recruited_bam"

    ###########################################################################
    # Convert back to paired FASTQs
    ###########################################################################

    samtools view \
        --threads 4 \
        -u \
        -f 0x1 \
        -F 0xF0C \
        "$subsampled_bam" |
    samtools collate \
        --threads 4 \
        -Ou \
        - |
    samtools fastq \
        --threads 4 \
        -n \
        -1 "$r1" \
        -2 "$r2" \
        -0 "$category0" \
        -s "$singleton" \
        -

    r1_reads=$(gzip -cd "$r1" | awk 'END { print int(NR / 4) }' )
    r2_reads=$(gzip -cd "$r2" | awk 'END { print int(NR / 4) }' )

    printf \
        "[INFO] %s: %d paired templates and %d singletons written\n" \
        "$sample" \
        "$r1_reads" \
        >&2

    rm -f \
        "$recruited_r1" \
        "$recruited_r2" \
        "$recruited_singleton" \
        "$recruited_bam" \
        "${recruited_bam}.bai" \
        "$subsampled_bam" \
        "${subsampled_bam}.bai"
done

# Create samplesheet

declare -A SAMPLE_IDS=(
    ["EM6"]="EM6"
    ["EM3"]="EM3"
    ["F3"]="F3A"
    ["F2xM12-F1"]="F2xM12-F1"
)

declare -A POPULATIONS=(
    ["EM6"]="Pop1"
    ["EM3"]="Pop1"
    ["F3"]="Pop2"
    ["F2xM12-F1"]="Pop3"
)

SAMPLESHEET="${OUTPUT_DIR}/test_samplesheet.csv"

printf 'sample,pop,fwd,rev\n' > "$SAMPLESHEET"

for sample in "${SAMPLES[@]}"; do
    sample_id="${SAMPLE_IDS[$sample]}"
    r1="${OUTPUT_DIR}/${sample}_subset_R1.fastq.gz"
    r2="${OUTPUT_DIR}/${sample}_subset_R2.fastq.gz"

    if [[ -z "${sample_id:-}" ]]; then
        echo "ERROR: no sample ID defined for source sample: $sample" >&2
        exit 1
    fi

    if [[ ! -v "POPULATIONS[$sample]" ]]; then
        echo "ERROR: no population defined for sample: $sample" >&2
        exit 1
    fi

    for fastq in "$r1" "$r2"; do
        if [[ ! -s "$fastq" ]]; then
            echo "ERROR: FASTQ is missing or empty: $fastq" >&2
            exit 1
        fi
    done

    printf '%s,%s,%s,%s\n' \
        "$sample_id" \
        "${POPULATIONS[$sample]}" \
        "$r1" \
        "$r2" \
        >> "$SAMPLESHEET"
done


```

### Run test datasets


Run the Qfly test dataset using the test profile
```
module purge
export NXF_VER=26.07.0-edge
module load Java/17

# Local execution on a BASC node using Conda
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR="/group/pathogens/IAWS/Personal/Alexp/conda_cache"
nextflow run . -profile conda,debug,test -resume

# BASC SLURM execution using Conda
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR="/group/pathogens/IAWS/Personal/Alexp/conda_cache"
nextflow run . -profile debug,test -config conf/basc.config --slurm_account fruitfly -resume


# Test BMSB
module purge
export NXF_VER=26.07.0-edge
module load Java/17
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR="/group/pathogens/IAWS/Personal/Alexp/conda_cache"
nextflow run . \
    -profile debug \
    -config conf/basc.config \
    --slurm_account fruitfly \
    --samplesheet sample_sheet.csv \
    --ref_genome /group/referencedata/mspd-db/genomes/insect/halyomorpha_halys/GCF_000696795.3_Hhal_1.1_genomic.fna \
    --mito_contig NC_013272.1 \
    -resume \
    -with-trace trace.txt scratch=false

```


# View cpu utilisation of running job

```
# List all jobs
squeue     --user "$USER"     --format='%i|%T|%j|%R' | grep 'MAP_TO_GENOME'

# Slect job_id
job_id=21074658

node=$(squeue -h -j "$job_id" -o '%N')

echo "Monitoring job ${job_id} on ${node}"

ssh "$node" bash -s -- "$job_id" <<'REMOTE'
job_id="$1"

pids=$(
    scontrol listpids "$job_id" 2>/dev/null |
        awk 'NR > 1 && $1 ~ /^[0-9]+$/ {print $1}' |
        paste -sd, -
)

if [[ -z "$pids" ]]; then
    echo "No processes found for Slurm job $job_id" >&2
    exit 0
fi

{
    printf '%-20s %10s %12s\n' \
        'COMMAND' \
        'CPU' \
        'RSS_GiB'

    ps -p "$pids" -o comm=,%cpu=,rss= |
        awk '
            {
                cpu[$1] += $2
                rss[$1] += $3
            }

            END {
                for (command in cpu) {
                    printf "%-20s %9.1f%% %12.2f\n",
                        command,
                        cpu[command],
                        rss[command] / 1024 / 1024
                }
            }
        ' |
        sort -k2,2nr
}
REMOTE
```