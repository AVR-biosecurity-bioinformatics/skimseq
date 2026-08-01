# skimseq
`skimseq` is a Nextflow-based bioinformatics pipeline for analysing low- and variable-coverage whole-genome sequencing data, with a particular focus on population genomics and large cohort studies. This pipeline is being developed by a team at [Agriculture Victoria Science and Technology](https://agriculture.vic.gov.au/)

> [!WARNING]
> This pipeline is currently under active development. The code, configuration and outputs may change without notice, and there is no guarantee that the current version is stable or suitable for production analyses.

# Background

Low-coverage whole-genome sequencing (often referred to as *skim sequencing* or *genome skimming*) has become a popular approach for population genomics because it enables large numbers of individuals to be sequenced at a fraction of the cost of traditional high-coverage whole-genome sequencing.

Many studies also incorporate publicly available sequencing data or samples generated across different projects, sequencing centres, library preparation methods, and platforms. As a result, coverage is often highly variable both within and between datasets. Even within a single study, differences in library quality, sequencing yield, and missing data can result in substantial variation in sequencing depth among samples.

A key challenge is that many variant-calling, filtering, and quality-control workflows were developed for relatively uniform, high-coverage datasets (typically 20-30× coverage). Applying the same assumptions and thresholds to low- and variable-coverage data can introduce bias, disproportionately exclude poorly covered samples, or retain low-quality variants in well-covered samples.

`skimseq` is designed specifically for low- and variable-coverage whole-genome sequencing data, providing an end-to-end workflow for read processing, alignment, variant discovery, quality filtering, and downstream population genomic analyses.

# Implementation

`skimseq` is a nextlow workflow that primarily relies on widely used community tools for sequence alignment and variant calling, with variant discovery performed using a `bcftools mpileup` and `bcftools call`-based workflow. These standard approaches are complemented by additional quality control, filtering, and masking procedures designed specifically for low- and variable-coverage datasets. Particular emphasis is placed on minimising biases associated with missing data, uneven coverage, repetitive regions, and problematic genomic sites.

`skimseq` broadly follows the nf-core template and conventions but intentionally diverges in several areas to prioritise computational efficiency, scalability, and reduced storage requirements. Many Nextflow workflows generate large numbers of intermediate files and repeatedly write data to disk. In contrast, `skimseq` aims to minimise I/O and temporary storage by streaming data between tools, combining logically related operations into single processes where appropriate, and avoiding unnecessary intermediate outputs.

As a result, some implementation details differ from typical nf-core recommendations. These design choices are intentional and are aimed at improving performance on shared high-performance computing (HPC) systems and large sequencing datasets while maintaining reproducibility and portability.

## Requirements

The pipeline requires:

- Linux, macOS or Windows using WSL 2
- Bash
- Java 17
- Nextflow
- one supported software backend:
  - Conda
  - Docker
  - Apptainer
  - Singularity

Conda dependencies are defined separately for each module using an `environment.yml` file. Container profiles use [Wave](https://docs.seqera.io/nextflow/wave) to provision containers from these Conda environment specifications.

## Installation

Clone the repository:

```
git clone https://github.com/AVR-biosecurity-bioinformatics/skimseq.git
cd skimseq
```

Install Nextflow if it is not already available
```
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow

mkdir -p "$HOME/.local/bin"
mv nextflow "$HOME/.local/bin/"

export PATH="$HOME/.local/bin:$PATH"
```

# Configuration profiles
Profiles are composable and should generally be supplied in the following order:
`execution backend, software backend, optional development profiles`

## Execution profiles:
- `local`: execute processes directly on the current machine
- `basc_slurm`: submit processes to BASC using SLURM
- **TODO** Add support for other schedulers

## Software profiles:
- `conda`: use module-level Conda environments
- `docker`: use Docker with Wave-provisioned containers
- `singularity`: use Singularity with Wave-provisioned containers **UNTESTED**
- `apptainer`: use Apptainer with Wave-provisioned containers **UNTESTED**
- `podman`: use Podman with Wave-provisioned containers **UNTESTED**

## Development profiles:
- test: use the bundled minimal Queensland fruit fly test dataset and reduced resources
- debug: retain intermediate process outputs and R session data
The test profile should be specified last so that its reduced resource limits override production settings.

# Usage examples:

## Running on AgVic BASC HPC using SLURM

The currently supported approach for running on AgVic BASC HPC uses the basc_conda software profile, and requires setting --slurm-account

```
# Load modules for launching job
module purge
module load Java/17
module load Miniconda3/24.7.1-0

# Export variables for nextflow:
export NXF_VER=26.07.0-edge
export NXF_CONDA_CACHEDIR="YOUR_CACHEDIR"

# Launch the job
nextflow run . \
    -profile basc \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

## Running a test run locally on a BASC node
For small tests, use the local executor instead of SLURM.
```
# Load modules for launching job
module purge
module load Java/17
module load Miniconda3/24.7.1-0

# Export variables for nextflow:
export NXF_VER=26.07.0-edge
export NXF_CONDA_CACHEDIR="YOUR_CACHEDIR"

# Launch the job
nextflow run . \
    -profile local,conda,debug,test \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

## Running locally with WSL and Conda
This assumes you already have WSL and a linux distribution installed (i.e. Ubuntu-24.04), as well as Nextflow and Miniconda

```
nextflow run . \
    -profile local,conda \
    -resume
```

# Setting pipeline parameters
Create a YAML parameter file containing the analysis inputs and settings, or edit the default parameters file found in `/assets/default_params.yml`
```
samplesheet: "/path/to/samplesheet.csv"
ref_genome: "/path/to/reference.fa"
mito_contig: "MT"

variant_caller: "bcftools"
calling_model: "cohort"
ploidy: 2

output_cram: true
output_gvcf: true
```
Supply it using `params-file`
```
# Launch the job
nextflow run . \
    -profile basc_slurm,basc_conda \
    --slurm_account YOUR_ACCOUNT \
    -params-file params.yml \
    -resume
```