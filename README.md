# skimseq
A nextflow-based bioinformatics pipeline for analysing low or variable coverage whole genome sequencing data

This pipeline is being developed by a team at [Agriculture Victoria Science and Technology](https://agriculture.vic.gov.au/), nd is currently intended primarily for internal use on Agriculture Victoria's BASC computing cluster.

> [!WARNING]
> This pipeline is currently underr active development. The code, configuration and outputs may change without notice, and there is no guarantee that the current version is stable or suitable for production analyses.

## Requirements

The pipeline requires:

- Linux, macOS or Windows using WSL 2
- Bash
- Java 17
- Nextflow
- one supported software backend:
  - Conda
  - centrally installed BASC modules
  - Docker
  - Apptainer
  - Singularity
  - Shifter
  - Podman
  - Charliecloud

Conda dependencies are defined separately for each module using an `environment.yml` file. Container profiles can use Wave to provision containers from the same Conda environment specifications.

## Installation

Clone the repository:

```
git clone https://github.com/AVR-biosecurity-bioinformatics/skimseq.git
cd skimseq
```

Install Nextflow if it is not alraedy available
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
- local: execute processes directly on the current machine
- basc_slurm: submit processes to BASC using SLURM

## Software profiles:
- basc_modules: use centrally installed software modules **FOR BASC SYSTEM ONLY**
- basc_conda: use module-level Conda environments **FOR BASC SYSTEM ONLY**
- conda: use module-level Conda environments
- docker: use Docker with Wave-provisioned containers **UNTESTED**
- apptainer: use Apptainer with Wave-provisioned containers **UNTESTED**
- singularity: use Singularity with Wave-provisioned containers **UNTESTED**
- podman: use Podman with Wave-provisioned containers **UNTESTED**
- shifter: use Shifter with Wave-provisioned containers **UNTESTED**
- charliecloud: use Charliecloud containers with Wave-provisioned containers **UNTESTED**

## Development profiles:
- test: use the bundled minimal Queensland fruit fly test dataset and reduced resources
- debug: retain intermediate process outputs and R session data

The test profile should be specified last so that its reduced resource limits override production settings.

# Running on AgVic BASC HPC using SLURM

The currently supported approach for running on BASC uses the conda software profile, and requires setting --slurm-account

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
    -profile basc_slurm,basc_conda \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

# Running a test run locally on a BASC node
For small tests within an interactive BASC allocation, use the local executor instead of SLURM.
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
    -profile local,basc_conda,debug,test \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

# Running locally with WSL and Conda
This assumes you already have WSL and a linux distribution installed (i.e. Ubuntu-24.04), as well as Nextflow and Miniconda

```
nextflow run . \
    -profile local,conda \
    -resume
```

# Setting pipeline parameters
Create a YAML parameter file containing the analysis inputs and settings:
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