#!/bin/bash
#SBATCH --job-name=vcfrelate         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB 
#SBATCH --time=6:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
#SBATCH --partition=shortrun
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                               SLURM Job handling                             -
#--------------------------------------------------------------------------------

# Only run job if it is submitted as a SLURM array
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=vcfrelate_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
vcf=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R Reference Genome ] [ -O output directory ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":V:S:O:t" options; do       
  # use silent error checking;
  case "${options}" in
    V)                             
      vcf=${OPTARG}
	  # Test if exists
	  if [ ! -f "$vcf" ] ; then  
        echo "Error: -V ${vcf} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo vcf=${vcf}	  
    ;;
    S)
	  sitelist=${OPTARG}
	  # Test if not empty
	  if [[ $sitelist ]]; then
		  # test if exists
		  if [ ! -f "$sitelist" ] ; then  
			echo "Error: -S ${sitelist} doesnt exist"
			exit_abnormal
			exit 1
		  fi
		echo sitelist=${sitelist}
	  fi
    ;;
    O)
      outdir=${OPTARG}
      echo outdir=${outdir}
    ;;
    t)
	  echo "-t has been set, subsampling all fastqs for testing purposes" >&2
	  test_run='true'
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

# Read input params
bamlist=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo bamlist=${bamlist}

# Get sample name from beagle file
Sample=$(basename ${vcf} .vcf.gz)
echo Sample=${Sample}

#set outdir file name
outname=$(echo $Sample $(echo $(basename ${bamlist}) | sed 's/bamlist_//g' | sed 's/.txt//g') $(basename ${sitelist} | sed 's/.sites.*$//g' ) | sed 's/ /-/g' )
echo outname=${outname}

# Make directories for outdirs
mkdir -p ${outdir}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

#--------------------------------------------------------------------------------
#-                              Subset indviduals and sites                       -
#--------------------------------------------------------------------------------
# Load BCFtools for filteering
module purge
module load BCFtools/1.21-GCC-13.3.0

# Get samples common to both sample_data and VCF file
comm -12  <(bcftools query -l ${vcf} | sort) <(cat ${bamlist} | sort) > samples_to_keep
echo $(cat samples_to_keep | wc -l ) samples

if [ -n "$sitelist" ]; then
	echo "Subsetting to only those sites within sitelist: ${sitelist}"
	if [[ "${sitelist}" == *.gz ]]; then
		zcat ${sitelist} | tr ':' '\t' > sites_to_keep
	else
		cat ${sitelist} | tr ':' '\t' > sites_to_keep
	fi
else 
	echo "no sitelist provided, using all sites"
	bcftools query -f '%CHROM\t%POS\n' ${vcf} > sites_to_keep
fi

# Print number of sites
echo $(cat sites_to_keep | wc -l) sites

# Subset  VCF to target samples and sites
bcftools view ${vcf} -S samples_to_keep -T sites_to_keep --threads ${SLURM_CPUS_PER_TASK} -Ou -o ${Sample}.vcf

#--------------------------------------------------------------------------------
#-                                    Run VCFtools                              -
#--------------------------------------------------------------------------------
#Load Modules
module purge
module load VCFtools/0.1.16-GCC-13.3.0

vcftools --vcf ${Sample}.vcf --out ${outname} --relatedness 

echo run ${outname} complete

# Copy files back to drive
cp ${outname}.relatedness ${outdir}/${outname}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 