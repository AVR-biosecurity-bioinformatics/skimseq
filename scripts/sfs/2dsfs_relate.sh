#!/bin/bash
#SBATCH --job-name=2dsfs         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB 
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

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=2dsfs_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi


#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
ReferenceGenome=""
AncestralGenome=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -O output directory ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts "R:A:O:" options; do       
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
	  # Test if not empty
	  if [[ $AncestralGenome ]]; then
		  # Test if exists
		  if [ ! -f "$AncestralGenome" ] ; then  
			echo "Error: -A ${AncestralGenome} doesnt exist"
			exit_abnormal
			exit 1
		  fi
	  	  echo AncestralGenome=${AncestralGenome}	  
	  fi
    ;;
    O)
      outdir=${OPTARG}
      echo outdir=${outdir}
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
# Input information 
comparisons=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo comparisons=${comparisons}
saf1=$(echo ${comparisons} | awk -F ':' '{print $1}' )
saf2=$(echo ${comparisons} | awk -F ':' '{print $2}' )

pop1=$(basename ${saf1} |  sed -r 's/.saf.idx//' )
pop2=$(basename ${saf2} |  sed -r 's/.saf.idx//' )

# Make outname
outname=$(echo ${pop1} ${pop2} | sed 's/ /-/g' )

# Make directories for outputs
mkdir -p ${outdir}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Check if ancestral genome is provided - if not then create folded sfs
if [[ $AncestralGenome ]]; then
	echo 'Ancestral genome provided, making unfolded sfs'

else
	echo 'No ancestral genome provided, making folded sfs'
	folded='true'
fi

## Copy saf files to temp
cp ${saf1} .
cp ${saf2} .
cp $( echo ${saf1} |  sed -r 's/.saf.idx/.saf.gz/' ) .
cp $( echo ${saf2} |  sed -r 's/.saf.idx/.saf.gz/' ) .
cp $( echo ${saf1} |  sed -r 's/.saf.idx/.saf.pos.gz/' ) .
cp $( echo ${saf2} |  sed -r 's/.saf.idx/.saf.pos.gz/' ) .

# Load modules 
module purge
module load HTSlib/1.21-GCC-13.3.0
module load parallel/20240722-GCCcore-13.3.0
module load Rust/1.78.0-GCCcore-13.3.0

#--------------------------------------------------------------------------------
#-                                   2D SFS & FST	                            -
#--------------------------------------------------------------------------------
echo "Calculating 2D SFS for ${pop1} and ${pop2}"

mkdir results

# Calculate 2D SFS from SAF files
# Use streaming mode to keep RAM low 
/home/ap0y/.cargo/bin/winsfs -v -t ${SLURM_CPUS_PER_TASK} ${pop1}.saf.idx ${pop2}.saf.idx > results/${outname}.sfs


if [[ "$folded" = true ]]; then
	/home/ap0y/.cargo/bin/winsfs view --fold results/${outname}.sfs > folded.sfs
	mv folded.sfs results/${outname}.sfs
fi

# Check the dimensions of the 2dsfs
dim1=$(grep "^#SHAPE=" results/${outname}.sfs | sed -E 's/^#SHAPE=<([0-9]+)\/([0-9]+)>/\1/')
dim2=$(grep "^#SHAPE=" results/${outname}.sfs | sed -E 's/^#SHAPE=<([0-9]+)\/([0-9]+)>/\2/')

echo "SFS dimension1: $dim1"
echo "SFS dimension2: $dim2"

if [[ "$dim1" -eq 3 && "$dim2" -eq 3 ]]; then
	echo "SFS is a 3x3 individual comparison, running individual stats" 
    # if 3x3 (2 individuals compared) - additionally calculate kinship statistics
	cat results/${outname}.sfs | /home/ap0y/.cargo/bin/sfs stat --statistics f2,fst,pi-xy,s,sum,king,r0,r1 > global_estimate.txt
else	
	echo 'both dimensions must equal 3'
	#exit
fi

# Create bootstrapped SFS - 100 reps
/home/ap0y/.cargo/bin/winsfs split -v -t ${SLURM_CPUS_PER_TASK} --sfs results/${outname}.sfs -S 100 ${pop1}.saf.idx ${pop2}.saf.idx  > results/${outname}_bootstrap.sfs
awk '/^#SHAPE/ {file = "boot_" ++i ".sfs"} {print > file}' results/${outname}_bootstrap.sfs

# Calculate stats for each bootstrap, and then leave one out stats
echo 'block,f2,fst,pi-xy,s,sum,king,r0,r1' > results/${outname}_block_estimates.txt
echo 'block,f2,fst,pi-xy,s,sum,king,r0,r1' > results/${outname}_loo_estimates.txt
paste -d , <(echo 'global') <(cat global_estimate.txt) >> results/${outname}_block_estimates.txt
paste -d , <(echo 'global') <(cat global_estimate.txt) >> results/${outname}_loo_estimates.txt

ls boot_*.sfs | while read -r file; do
	block=$(echo $file | sed 's/.sfs//g')
	echo calculating bootstrap estimate for block "$block"
	
	# if ancestral genome wasnt provided, fold the bootstrapped sfs before calculations
	if [[ "$folded" == "true" ]]; then
        cat "$file" | /home/ap0y/.cargo/bin/winsfs view --fold > tmp.sfs
	else
        cp "$file" tmp.sfs
	fi
	
	# Run different stats depending on dimensions
	if [[ "$dim1" -eq 3 && "$dim2" -eq 3 ]]; then
        # if 3x3 (2 individuals compared) - additionally calculate kinship statistics
        cat tmp.sfs | /home/ap0y/.cargo/bin/sfs stat --statistics f2,fst,pi-xy,s,sum,king,r0,r1 > block_estimate.txt
	else	
		echo 'both dimensions must equal 3'
		#exit
	fi
	
	paste -d , <(echo $block) <(cat block_estimate.txt) >> results/${outname}_block_estimates.txt
	
	# Subtract block estimates from global estimates to get leave-one-out-blocks
	awk 'NR==FNR {split($0, block, ","); next} {for (i=1; i<=NF; i++) printf "%.6f%s", $i - block[i], (i<NF ? "," : "\n")}' block_estimate.txt FS="," OFS="," global_estimate.txt > loo.txt
	paste -d , <(echo $block) <(cat loo.txt) >> results/${outname}_loo_estimates.txt
	
	# Remove temporary file
	rm -f tmp.sfs
    rm -f block_estimate.txt
done


# Copy files back to drive
cp -r results/* ${outdir}/${outname}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
