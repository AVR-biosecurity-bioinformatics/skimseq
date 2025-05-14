#!/bin/bash
#SBATCH --job-name=eems         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=40GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

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

Index=eems_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

# Read input params
distmat=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo distmat=${distmat}

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -D distance matrix ] [ -S site list ] [ -H habitat coordinates ] [ -O output directory ]" 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
# Do i need to add manifest
while getopts "C:H:O:" options; do
  # use silent error checking;
  case "${options}" in
    C)                             
      coord=${OPTARG}
	  # Test if exists
	  if [ ! -f "$coord" ] ; then  
        echo "Error: -C ${coord} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo coord=${coord}	  
      ;;
    #S)
	#  sitelist=${OPTARG}
	#  # test if exists
	#  if [ ! -f "$sitelist" ] ; then  
	#	echo "Error: -S ${sitelist} doesnt exist"
	#	exit_abnormal
	#	exit 1
	#  fi
    #  echo sitelist=${sitelist}
    #;;
    H)                             
      habitat=${OPTARG}
	  # Test if exists
	  if [ ! -f "$habitat" ] ; then  
        echo "Error: -H ${habitat} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo habitat=${habitat}	  
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

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------

# Get sample name from distmat
Sample=$(basename ${coord} | sed 's/.tsv//g' | sed 's/.txt//g')
echo Sample=${Sample}

# sitelist
sitelist=$(echo $distmat | sed 's/.ibsMat.*$/.sitecounts.txt/g')

#set output file name
outname=$(echo $Sample $(echo  $(basename ${distmat} | sed 's/.ibsMat.*$//g' )) | sed 's/ /-/g' )

# Get the sample names for the input dist matrix
input_names=$(find $(dirname $distmat ) | grep _bams.txt)

# Make directories for outputs
mkdir -p ${output}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

#Load Modules
module purge
module load eems/20231029-GCC-13.3.0
module load R/4.4.2-gfbf-2024a

#--------------------------------------------------------------------------------
#-                                Rearrange inputs                              -
#--------------------------------------------------------------------------------

# Create R script 
echo -e "input_args <- commandArgs(TRUE)
library(tidyverse)

# Process input arguments
distmat_file <- input_args[1]
sample_names <- readLines(input_args[2])%>%
  str_remove('.bam')
coords <- read_tsv(input_args[3], col_names=c('long', 'lat', 'sample_id'))
habitat <- read_tsv(input_args[4], col_names=c('long', 'lat', 'sample_id'))

# Check coords
if(any(duplicated(coords))){
  dups <- coords\$sample_id[duplicated(coords\$sample_id)]
  stop(paste0('Duplicated sample IDs in coords file: ',paste(dups, collapse=', ')))
}

# Subset coords file to just those in habitat boundariesz
coords2 <- coords %>%
  dplyr::filter(between(lat, min(habitat\$lat), max(habitat\$lat))) %>%
  dplyr::filter(between(long, min(habitat\$long), max(habitat\$long)))

# Get list of samples that are in both input matrix and coords file
keeplist <- intersect(sample_names,coords2\$sample_id)

#Subset coords2 file to just those samples
coords2 <- coords2 %>% filter(sample_id %in% keeplist) %>%
  arrange(factor(sample_id, levels = keeplist)) %>%
  dplyr::select(-sample_id)
write_tsv(coords2, 'eems_input.coord', col_names = FALSE)

# Subset distmat
distmat <- read_tsv(distmat_file, col_names = sample_names) %>%
  select_if(function(x) !(all(is.na(x)) | all(x=='')))%>%
  as.data.frame() %>%
    magrittr::set_rownames(sample_names)

#print(distmat)
# subset distmat to just include those samples
distmat2 <- distmat[keeplist,keeplist]
write_tsv(distmat2, 'eems_input.diffs', col_names = FALSE)
message('file subset complete')
" > reformat_eems_inputs.R

Rscript --vanilla reformat_eems_inputs.R ${distmat} ${input_names} ${coord} ${habitat}

cp ${habitat} eems_input.outer

#--------------------------------------------------------------------------------
#-                                    Run EEMS                                  -
#--------------------------------------------------------------------------------

# Create parameters file for eems
nind=$(cat eems_input.diffs | wc -l)
nsites=$(cat ${sitelist} | wc -l)
echo $nind individuals and $nsites sites

# Check dimensions of distance matrix
ncol=$(awk '{print NF}' eems_input.diffs | sort -nu | tail -n 1)
nrow=$(cat eems_input.diffs | wc -l)
echo $ncol columns and $nrow rows in distance matrix

echo -e "datapath = eems_input
mcmcpath = eems_output
nIndiv = ${nind}
nSites = ${nsites}
nDemes = 200
diploid = true
numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999" > eems_params.txt

# Run eems using input params
runeems_snps --params eems_params.txt

# Copy files back to drive
cp eems_output/* ${output}/${outname}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
