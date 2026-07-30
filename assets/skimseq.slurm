#!/bin/bash

# The name of the job:
#SBATCH --job-name="skimseq"

# The project ID which this job should run under:
#SBATCH --account=

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# The total amount of memory in megabytes in the job:
#SBATCH --mem-per-cpu=8GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

#SBATCH --mail-user=

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=14-00:0:00

# Output errors and logs into same file
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

### run job
module load Java/17
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR=""
export NXF_VER=26.07.0-edge

nextflow run . -profile basc_slurm,basc_conda \
   -ansi-log false \
   --samplesheet sample_sheet.csv -c qfly.config -resume -with-trace trace.txt 
