#!/bin/bash
#SBATCH --job-name=locator  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB 
#SBATCH --time=6:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly2
#SBATCH --partition=shortrun,batch
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out


# ------------------------------------------------------------------------------
# Usage (array):
#   N=$(awk 'BEGIN{n=0} $0!~/^(#|track|browser)/ && NF>=3 {n++} END{print n}' regions.bed)
#   sbatch --array=1-${N}%20 run_locator_array.slurm \
#       -V input.vcf.gz -S sample_data.txt -O results_dir -B regions.bed \
#       [-T sites.txt]
# ------------------------------------------------------------------------------

USAGE="Usage:
  sbatch --array=1-N $(basename "$0") \\
    -V input.vcf.gz -S sample_data.txt -O outdir -B regions.bed \\
    [--min-mac 2] [--impute-missing] [-T sitelist]

Required:
  -V, --vcf            Input VCF(.gz or .bcf)
  -S, --sample-data    Sample data file (x, y, sampleID)
  -O, --outdir         Output directory
  -B, --regions-bed    BED file of regions (CHROM<TAB>START<TAB>END), 0-based half-open. Array index selects the Nth region.

Parameters:
  --min-mac N          Minimum minor allele count for filtering before training (default: 2)

Optional:
  --impute-missing     Impute missing genotypes (locator flag; default: off)
  -T, --sites          Site list file (CHROM:POS or CHROM<TAB>POS); will be intersected with region
"


die_usage() {
    [[ $# -gt 0 ]] && echo "Error: $1" >&2
    printf '%s\n' "$USAGE" >&2
    exit 1
}

# Default params
VCF=""
SAMPLEDATA=""
OUTDIR=""
SITELIST=""
REGIONS_BED=""
MIN_MAC=2
IMPUTE_MISSING=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -V|--vcf)         [[ $# -ge 2 ]] || die_usage "Missing value for $1"; VCF="$2"; shift 2 ;;
        -S|--sample-data) [[ $# -ge 2 ]] || die_usage "Missing value for $1"; SAMPLEDATA="$2"; shift 2 ;;
        -O|--outdir)      [[ $# -ge 2 ]] || die_usage "Missing value for $1"; OUTDIR="$2"; shift 2 ;;
        -T|--sites)       [[ $# -ge 2 ]] || die_usage "Missing value for $1"; SITELIST="$2"; shift 2 ;;
        -B|--regions-bed) [[ $# -ge 2 ]] || die_usage "Missing value for $1"; REGIONS_BED="$2"; shift 2 ;;
        --min-mac)        [[ $# -ge 2 ]] || die_usage "Missing value for $1"; MIN_MAC="$2"; shift 2 ;;
        --impute-missing) IMPUTE_MISSING=1; shift ;;

        -h|--help)        die_usage; exit 0 ;;
        *)                die_usage "Unknown argument: $1" ;;
    esac
done

# Validate parameters
[[ -n "$VCF" ]] || die_usage "Missing -V/--vcf"
[[ -n "$SAMPLEDATA" ]] || die_usage "Missing -S/--sample-data"
[[ -n "$OUTDIR" ]] || die_usage "Missing -O/--outdir"
[[ -n "$REGIONS_BED" ]] || die_usage "Missing -B/--regions-bed"


if [[ ! "${MIN_MAC}" =~ ^[0-9]+$ ]]; then
    die_usage "--min-mac must be a non-negative integer, got: ${MIN_MAC}"
fi
if (( MIN_MAC < 0 )); then
    die_usage "--min-mac must be >= 0, got: ${MIN_MAC}"
fi


# Require array context
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a SLURM array (SLURM_ARRAY_TASK_ID not set)}"

mkdir -p "${OUTDIR}"


# Convert relative paths to absolute
VCF=$(readlink -f "$VCF")
SAMPLEDATA=$(readlink -f "$SAMPLEDATA")
OUTDIR=$(readlink -f "$OUTDIR")
REGIONS_BED=$(readlink -f "$REGIONS_BED")

if [[ -n "$SITELIST" ]]; then
    [[ -f "$SITELIST" ]] || die_usage "Site list not found: $SITELIST"
    SITELIST=$(readlink -f "$SITELIST")
fi


# ------------------------------------------------------------------------------
# Select the array task's region from BED (ignores track/browser/# lines)
# ------------------------------------------------------------------------------
REG_LINE=$(awk -v n="${SLURM_ARRAY_TASK_ID}" '
    BEGIN{FS="\t"; i=0}
    $0 ~ /^(#|track|browser)/ {next}
    NF < 3 {next}
    {i++; if (i==n) {print $1, $2, $3; exit}}
' "${REGIONS_BED}")

[[ -n "${REG_LINE}" ]] || { echo "No region found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

CHR=$(echo "${REG_LINE}" | awk '{print $1}')
START0=$(echo "${REG_LINE}" | awk '{print $2}')
END0=$(echo "${REG_LINE}" | awk '{print $3}')

# BED (0-based, half-open) -> bcftools (1-based, inclusive)
FROM=$((START0 + 1))
TO=$((END0))
REGION="${CHR}:${FROM}-${TO}"

# Region tag for naming
REGION_TAG="${CHR}_${START0}_${END0}"


# ------------------------------------------------------------------------------
# RUNNAME includes region so outputs don't collide
# ------------------------------------------------------------------------------
BASENAME=$(basename "${VCF}")
BASENAME=${BASENAME%.vcf.gz}; BASENAME=${BASENAME%.vcf}; BASENAME=${BASENAME%.bcf}

RUNNAME="${BASENAME}_${REGION_TAG}"
if [[ -n "${SITELIST}" ]]; then
    RUNNAME="${RUNNAME}-$(basename "${SITELIST}" | sed 's/\.[^.]*$//')"
fi

# ------------------------------------------------------------------------------
# Workdir (include array task to avoid collisions)
# ------------------------------------------------------------------------------
WORKDIR="${TMPDIR:-/tmp}/${USER}/locator_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "WORKDIR=${WORKDIR}"
echo "RUNNAME=${RUNNAME}"
echo "VCF=${VCF}"
echo "SAMPLEDATA=${SAMPLEDATA}"
echo "OUTDIR=${OUTDIR}"
echo "SITELIST=${SITELIST:-}"
echo "REGIONS_BED=${REGIONS_BED}"
echo "REGION=${REGION}"


# ------------------------------------------------------------------------------
# Load software
# ------------------------------------------------------------------------------
module purge
module load BCFtools/1.23.1-GCC-13.3.0

# conda init for non-interactive shell
source /group/pathogens/IAWS/Personal/Alexp/miniconda3/etc/profile.d/conda.sh
conda activate locator

echo "Python: $(which python)"
echo "Locator: $(which locator)"
echo "bcftools: $(which bcftools)"

# ------------------------------------------------------------------------------
# Determine overlapping samples between VCF and sample_data
# ------------------------------------------------------------------------------
bcftools query -l "${VCF}" | sort > vcf.samples.txt
awk 'NR>1 {gsub(/"/, "", $3); print $3}' "${SAMPLEDATA}" | sort > sample_data.samples.txt
comm -12 vcf.samples.txt sample_data.samples.txt > samples_to_keep.txt

NSAMPLES=$(wc -l < samples_to_keep.txt)
echo "${NSAMPLES} overlapping samples found"

if [[ "${NSAMPLES}" -eq 0 ]]; then
    echo "No overlapping samples between VCF and sample_data"
    exit 1
fi

# filter sample_data to overlapping samples
awk 'BEGIN{FS=OFS="\t"}
     FNR==NR {keep[$1]=1; next}   # first file: samples_to_keep.txt
     FNR==1  {print; next}        # second file: print header of sample_data
     ($3 in keep) {print}' \
     samples_to_keep.txt "${SAMPLEDATA}" > "${RUNNAME}_sample_data.txt"


# ------------------------------------------------------------------------------
# Prepare sites file if requested
# ------------------------------------------------------------------------------
BCFTOOLS_SITE_ARG=()
if [[ -n "${SITELIST}" ]]; then
    [[ -f "${SITELIST}" ]] || { echo "Sitelist not found: ${SITELIST}" >&2; exit 1; }

    # Decompress/copy to raw text
    if [[ "${SITELIST}" == *.gz ]]; then
        zcat "${SITELIST}" > raw_sites.txt
    else
        cp "${SITELIST}" raw_sites.txt
    fi

    # Convert to BED (0-based, half-open): chrom  start=pos-1  end=pos
    # Accepts either "CHROM:POS" or "CHROM<TAB>POS"
    awk 'BEGIN{FS=OFS="\t"}
         NF==0 {next}
         {
             chrom=""; pos=""
             if (index($0,":")>0) {
                 split($0,a,":"); chrom=a[1]; pos=a[2]
             } else {
                 chrom=$1; pos=$2
             }
             # basic numeric check
             if (pos !~ /^[0-9]+$/) next

             start=pos-1
             if (start < 0) start=0
             end=pos
             print chrom, start, end
         }' raw_sites.txt > sites_to_keep.bed

    # Sort for tabix (required)
    LC_ALL=C sort -k1,1 -k2,2n sites_to_keep.bed > sites_to_keep.sorted.bed

    NSITES=$(wc -l < sites_to_keep.sorted.bed)
    echo "${NSITES} target sites provided"

    # bgzip + tabix index (BED preset)
    bgzip -c sites_to_keep.sorted.bed > sites_to_keep.bed.gz
    tabix -f -p bed sites_to_keep.bed.gz

    # Use -R with tabixed regions file
    BCFTOOLS_SITE_ARG=(-R sites_to_keep.bed.gz)
else
    echo "No sitelist provided; using all sites"
fi


# ------------------------------------------------------------------------------
# Region arg for this array task 
# ------------------------------------------------------------------------------
BCFTOOLS_REGION_ARG=(-r "${REGION}")

# ------------------------------------------------------------------------------
# Build minimal header matching subset samples + region + (optional) sites
# ------------------------------------------------------------------------------
bcftools view -h \
    -S samples_to_keep.txt \
    "${BCFTOOLS_REGION_ARG[@]}" \
    "${BCFTOOLS_SITE_ARG[@]}" \
    "${VCF}" \
| awk '
  BEGIN { kept_fileformat=0; kept_gt=0 }
  /^##fileformat=/ { if (!kept_fileformat) { print; kept_fileformat=1 } ; next }
  /^##contig=/     { print; next }
  /^##INFO=/       { print; next }
  /^##FORMAT=<ID=GT,/ { if (!kept_gt) { print; kept_gt=1 } ; next }
  /^#CHROM/        { print; next }
' > minimal.hdr

# ------------------------------------------------------------------------------
# Subset VCF to region + samples (+ optional sites), keep GT only, reheader, write
# ------------------------------------------------------------------------------
bcftools view -Ou \
    "${VCF}" \
    -S samples_to_keep.txt \
    "${BCFTOOLS_REGION_ARG[@]}" \
    "${BCFTOOLS_SITE_ARG[@]}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
| bcftools annotate -Ou \
    -x ID,QUAL,INFO,^FORMAT/GT \
| bcftools reheader -h minimal.hdr \
| bcftools view -Oz --threads "${SLURM_CPUS_PER_TASK}" \
    -o "${RUNNAME}.vcf.gz"

tabix -p vcf "${RUNNAME}.vcf.gz"

N_SAMPLES=$(bcftools query -l "${RUNNAME}.vcf.gz" | wc -l)
N_VARIANTS=$(bcftools index -n "${RUNNAME}.vcf.gz" 2>/dev/null || echo 0)

echo "${N_SAMPLES} samples after subsetting"
echo "${N_VARIANTS} variants after subsetting"

if (( N_VARIANTS == 0 )); then
    echo "ERROR: No variants remain after subsetting." >&2
    echo "       RUNNAME=${RUNNAME}" >&2
    [[ -n "${REGION:-}" ]] && echo "       REGION=${REGION}" >&2
    [[ -n "${SITELIST:-}" ]] && echo "       SITELIST=${SITELIST}" >&2
    echo "       This usually means the region/window contains no SNPs or filters removed everything." >&2
    exit 2
fi

# ------------------------------------------------------------------------------
# Convert to zarr
# ------------------------------------------------------------------------------
vcf_to_zarr --vcf "${RUNNAME}.vcf.gz" --zarr "${RUNNAME}.zarr" --overwrite

# ------------------------------------------------------------------------------
# Run locator
# ------------------------------------------------------------------------------

# Locator optional args

LOCATOR_ARGS=(--min_mac "${MIN_MAC}")
if (( IMPUTE_MISSING == 1 )); then
    LOCATOR_ARGS+=(--impute_missing)
fi


# Note, needs to be output in 'out' directory or error occurs
mkdir -p out

locator \
  --zarr "${RUNNAME}.zarr" \
  --sample_data "${RUNNAME}_sample_data.txt" \
  --out out/${RUNNAME} \
    "${LOCATOR_ARGS[@]}" \
  --disable_gpu \
  --keras_verbose 0 \
  --verbose
  
# ------------------------------------------------------------------------------
# Create custom output file
# ------------------------------------------------------------------------------
  
# Calculate missing data (per sample) + total records as a column
MISSING_TSV="${RUNNAME}.missing.tsv"

bcftools stats --threads "${SLURM_CPUS_PER_TASK}" -s - "${RUNNAME}.vcf.gz" \
| awk -v out="${MISSING_TSV}" 'BEGIN{OFS="\t"; FS="\t"}
    # capture PSC header and locate nMissing + sample columns dynamically
    $1=="#" && $2=="PSC" {
        for (i=3; i<=NF; i++) {
            if ($i=="sample")   sample_col=i
            if ($i=="nMissing") miss_col=i
        }
        next
    }

    $1=="PSC" {
        if(!printed_header){
            print "SAMPLE","NMISS" > out
            printed_header=1
        }
        s = (sample_col ? $(sample_col) : $3)
        m = (miss_col   ? $(miss_col)   : $14)
        print s, m > out
    }
'


# Find the predlocs file for this run
PREDLOCS_FILE=$(ls -1 "out/${RUNNAME}_predlocs.txt" 2>/dev/null | head -n 1)

PREDLOCS_TSV="tmp.tsv"

awk -F',' -v OFS='\t' '{
  # basic CSV->TSV (safe here because predlocs has no quoted commas)
  print $1, $2, $3
}' "${PREDLOCS_FILE}" > "${PREDLOCS_TSV}"


PREDLOCS_TSV_OUT="out/${RUNNAME}_predlocs.tsv"

awk -v FS='\t' -v OFS='\t' -v nsites="${N_VARIANTS}" '
  NR==FNR {
    # first file: missing TSV
    if (FNR==1) next
    miss[$1] = $2
    next
  }
  FNR==1 {
    # second file: predlocs TSV header
    # standardise to sampleID (optional)
    print $0, "n_genotyped", "n_sites"
    next
  }
  {
    id = $1
    m = (id in miss ? miss[id] : nsites)  # if no entry, assume fully missing
    g = nsites - m
    if (g < 0) g = 0
    print $0, g, nsites
  }
' "${MISSING_TSV}" "${PREDLOCS_TSV}" > "${PREDLOCS_TSV_OUT}"

rm "${PREDLOCS_FILE}" "${PREDLOCS_TSV}"

# ------------------------------------------------------------------------------
# Run VCF2Dis on the same chunk
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Copy outputs back
# ------------------------------------------------------------------------------
mkdir -p "${OUTDIR}/."
cp -r out/* "${OUTDIR}/."

echo "Completed successfully"

# Output useful job stats
/usr/local/bin/showJobStats.scr 