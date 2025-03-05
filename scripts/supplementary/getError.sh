#!/bin/bash
#SBATCH --job-name=angsdErr       
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                        Install sofware from github                           -
#--------------------------------------------------------------------------------

# Install latest development version of software into virtual environment
# module purge
# module load Python/3.8.2-GCCcore-9.3.0 
# module load GSL/2.7-GCC-11.2.0
# module load HTSlib/1.10.2-GCC-8.2.0-2.31.1
# virtualenv ~/angsd
# source ~/angsd/bin/activate
# cd angsd
# git clone https://github.com/angsd/angsd.git;
# cd angsd;make
# cd ..
# git clone https://github.com/fgvieira/ngsLD.git
# cd ngsLD;make
# cd ..
# git clone https://github.com/fgvieira/prune_graph.git
# cd prune_graph;cargo build --release

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------

# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

# This script calls genotype likelihoods in merged BAMs using ANGSD

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=fasta_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
ReferenceGenome=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R ReferenceGenome ] [ -t ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":R:A:I:O:" options; do       
  # use silent error checking;
  case "${options}" in
    R)                             
      ReferenceGenome=${OPTARG}
	  # Test if exists
	  if [ ! -f "$ReferenceGenome" ] ; then  
        echo "Error: -R ${ReferenceGenome} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo ReferenceGenome=${ReferenceGenome}	  
      ;;
    A)                             
      AncestralGenome=${OPTARG}
	  # Test if exists
	  if [ ! -f "$AncestralGenome" ] ; then  
        echo "Error: -A ${AncestralGenome} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo AncestralGenome=${AncestralGenome}	  
    ;;
    I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;
    O)
      output=${OPTARG}
      echo output=${output}
    ;;
	:) 
	# Exit If expected argument omitted
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal 
	  ;;
    *) 
	# Exit If unknown (any other) option
      exit_abnormal
      ;;
  esac
done
shift $((OPTIND -1))

# Exit if no ref genome supplied
if [[ "$ReferenceGenome" = "" ]]; then
	exit_abnormal
fi

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------
# Read sample from input index
Sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo Sample=${Sample}
echo SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}
echo SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS for input job
find $(/usr/bin/ls -d ${indir}) | grep -F ${Sample} | sort | uniq  > ${Sample}_tmp.txt
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_tmp.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

# Copy data files to temp and decompress
cp ${ReferenceGenome} .
sleep 20
cp ${ReferenceGenome}.fai .
cp ${AncestralGenome} .
sleep 20
cp ${AncestralGenome}.fai .
xargs -a ${Sample}_tmp.txt cp -ft .

#Load Modules
module purge
module load GSL/2.7-GCC-11.2.0
module load SAMtools/1.15-GCC-11.2.0
module load HTSlib/1.15-GCC-11.2.0
module load BCFtools/1.15-GCC-11.2.0
module load R/4.2.0-foss-2021b
#module load angsd/0.933-GCC-9.3.0

# Check BAMs
samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'

#--------------------------------------------------------------------------------
#-                                GetError                                      -
#--------------------------------------------------------------------------------

# launch angsd (local github version) 
# Using -rf to subset to target region
# Using only read filters - leave site filters for followup script
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 2 : realign around indels using method from original GATK paper
# -minQ 15: only bases with <0.03% prob of error # 
/home/ap0y/angsd/angsd/angsd -bam ${Sample}_bams.txt \
-ref $(basename ${ReferenceGenome}) -anc $(basename ${ReferenceGenome}) \
-minMapQ 20 -baq 2 -C 50 -minQ 15 \
-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
--ignore-RG 0 \
-doAncError 1 \
-nThreads ${SLURM_CPUS_PER_TASK} \
-out ${Sample} 

# Run the R script to calculate the errors
Rscript /home/ap0y/angsd/angsd/R/estError.R file=${Sample}.ancError indNames=${Sample}

# Copy files back to drive
cp -r ${Sample}.ancError ${output}/${Sample}.ancError
cp -r ${Sample}.ancErrorChr ${output}/${Sample}.ancErrorChr
cp -r errorEstOverChr.txt ${output}/${Sample}_errorEstOverChr.txt
cp -r errorEst.txt ${output}/${Sample}_errorEst.txt

# Output useful job stats
/usr/local/bin/showJobStats.scr 
