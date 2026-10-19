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

`skimseq` is a nextlow workflow that currently uses `bcftools mpileup` and `bcftools call`-based workflow. `skimseq` is developed for production-scale analyses and prioritises computational efficiency, scalability, and reduced storage requirements. In particular, `skimseq` aims to minimise I/O and temporary storage by streaming data between tools, combining logically related operations into single processes where appropriate, and avoiding unnecessary intermediate outputs. This approach diverges intentionally from standard nf-core pipelines.

## Requirements

The pipeline requires:

- Linux, macOS, or Windows using WSL 2
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

Install a local version of nextflow if it is not already available on your rsystem.
```
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow

mkdir -p "$HOME/.local/bin"
mv nextflow "$HOME/.local/bin/"

export PATH="$HOME/.local/bin:$PATH"
```

# Configuration profiles

Profiles are composable. In most cases, an institutional profile should be specified first to configure where tasks are executed, followed by a software environment profile specifying how process dependencies are provided. Optional testing and development profiles can be added last.

## Institutional profiles
Institutional profiles configure the execution environment for a specific compute system, including the scheduler/executor, queues, resource limits, filesystem settings, and other site-specific configuration.

Included profiles:
- `basc`: Agriculture Victoria BASC cluster (SLURM + Conda)
  - Requires `--slurm_account ACCOUNT`
- `gadi`: NCI Gadi cluster (PBS Pro + Singularity)
- config files for oher HPCs can be obtained from [nf-core/configs](https://github.com/nf-core/configs)

## Software environment profiles:

A software environment profile must be selected to provide the software dependencies required by each process.
- `conda`: create and use the module-level Conda environments directly
- `docker`: use Docker with Wave-provisioned containers
- `singularity`: use Singularity with Wave-provisioned containers
- `apptainer`: use Apptainer with Wave-provisioned containers
- `podman`: use Podman with Wave-provisioned containers

For the container profiles, Wave dynamically provisions the container images required by the workflow from the software dependencies defined by each module. Wave provides the container image, while the selected container runtime (Docker, Singularity, Apptainer, or Podman) is responsible for executing that image on the compute system.

The selected runtime must therefore be installed and supported on the target system. 

## Development profiles:
- `test`: use the bundled minimal Queensland fruit fly test dataset and reduced resources
  - The `test` profile should be specified last so that its reduced resource limits override production settings.

# Usage examples:

## Running on AgVic BASC HPC using SLURM

The currently supported approach for running on AgVic BASC HPC uses `conf/basc.config` institutional profile, and requires setting --slurm-account

```
# Load modules for launching job
module purge
module load Java/17
module load Miniconda3/24.7.1-0

# Launch the job
nextflow run . \
    -config conf/basc.config \
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

# Launch the job
nextflow run . \
    -profile conda,test \
    -resume
```

## Running locally with WSL and Docker Desktop
This assumes you already have WSL and a linux distribution installed (i.e. Ubuntu-24.04), as well as Nextflow and docker desktop

```
nextflow run . \
    -profile docker \
    -resume
```

## Running test dataset on NCI GADI HPC

**NOTE: THIS IS STILL TO BE TESTED**
```
# Clone the pipeline to your analysis directory
git clone https://github.com/AVR-biosecurity-bioinformatics/skimseq.git path/to/your/analysis/dir
cd path/to/your/analysis/dir

# Check your NCI gadi PROJECT and USER environmental variables are set (may need to set yourself?)
$PROJECT
$USER

# Check nextflow version
nextflow -v

# Run the pipeline - setting the work directory (temporary files) to the scratch directory
nextflow run . \
    -config conf/nci_gadi.config \
    -work-dir /scratch/${PROJECT}/${USER}/work \
    -profile test \
    -resume
```


# Setting pipeline parameters
Pipeline parameters can be supplied either in a YAML parameter file or directly on the command line. A parameter file is recommended for routine analyses because it provides a reproducible record of the inputs and settings used for a run.

## Create a parameters file
Create a YAML parameter file containing the analysis inputs and settings, or edit the example/default parameter file in assets/default_params.yml:
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

Supply the file using `-params-file`
```
# Launch the job
nextflow run . \
    -config conf/basc.config \
    -params-file params.yml \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

## Setting parameters on the command line

Pipeline parameters can also be provided directly on the command line using `--<parameter>`. For example, the equivalent inputs and calling settings could be specified as:

```
nextflow run . \
    -config conf/basc.config \
    --samplesheet /path/to/samplesheet.csv \
    --ref_genome /path/to/reference.fa \
    --mito_contig MT \
    --variant_caller bcftools \
    --calling_model cohort \
    --ploidy 2 \
    --output_cram true \
    --output_gvcf true \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

Command-line parameters can also be used together with -params-file, which is useful when most settings should come from a parameter file but one or two need to be changed for a particular run. Command line parameters are prioritised over those in the `params.yml` file

```
nextflow run . \
    -config conf/basc.config \
    -params-file params.yml \
    --output_cram false \
    --slurm_account YOUR_ACCOUNT \
    -resume
```

# Full list of available parameters

see `docs/parameters.md` for a full list of available pipeline parameters and their defaults.
This can also be seen in the `main.nf` file.
