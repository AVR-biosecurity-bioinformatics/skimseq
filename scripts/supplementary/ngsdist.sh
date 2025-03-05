#!/bin/bash
#SBATCH --job-name=ngsDist       
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=48:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------

# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

# This script uses NGSDist to estimate pairwise distances between individuals using genotype likelihoods

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=ngsdist_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi


#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
beagle=""

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
while getopts ":B:S:R:O:t" options; do       
  # use silent error checking;
  case "${options}" in
    B)                             
      beagle=${OPTARG}
	  # Test if exists
	  if [ ! -f "$beagle" ] ; then  
        echo "Error: -B ${beagle} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo beagle=${beagle}	  
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

# Get sample name from beagle file
Sample=$(basename ${beagle} .beagle.gz)
echo Sample=${Sample}

#set outdir file name
outname=$(echo $Sample $(echo $(basename ${bamlist}) | sed 's/bamlist_//g' | sed 's/.txt//g') $(basename ${sitelist} | sed 's/.sites.*$//g' ) | sed 's/ /-/g' )

# Make directories for outdirs
mkdir -p ${outdir}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Copy data files to temp and decompress
cp ${beagle} .
cp ${sitelist} .
cp ${bamlist} .

ls 
pigz -p ${SLURM_CPUS_PER_TASK} -d *.gz

#Load Modules
module purge
module load ngsDist/0f42786-GCC-9.3.0
module load FastME/2.1.6.1-GCCcore-8.2.0
module load RAxML-NG/0.9.0-GCCcore-8.2.0

#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------
# if -t is set, all fastqs will be subsampled to make the run fastQC
if [[ "$test_run" = true ]]; then
    # For quick testing - subsample input to just first chromosome and 1000 sitees
    echo subsampling to first chr/contig only
    chr1=$(cat ${Sample}.beagle | head -n 2 | tail -n 1 | awk -F '\t' '{print $1}' | awk -F '_' '{print $1}')
    cat ${Sample}.beagle | head -n 1  > ${Sample}.test.beagle
    cat ${Sample}.beagle | grep ${chr1} | head -1000 >> ${Sample}.test.beagle
    mv ${Sample}.test.beagle ${Sample}.beagle
fi

#--------------------------------------------------------------------------------
#-                           		 Prune samples                              -
#--------------------------------------------------------------------------------

# Subset beagle to only the samples within bamlist
cat ${Sample}.beagle | head -1 | tr "\t" "\n" > header

# Get file to search
printf  'marker\nallele1\nallele2\n' > search.txt
cat ${bamlist} | sed 's/.bam//g' >> search.txt
 
# Get the rows that arent in the bamlist
grep -vf search.txt header | sort | uniq > to_remove

#Check if any samples need removing
if grep -q '[^[:space:]]' to_remove; then
	echo "Subsetting beagle file to match bamlist"	
	# Get the index of the rows to remove and put in a variable
	index_rem=$(grep -nFf to_remove header | cut -d : -f 1 | tr "\n" ","| sed '$s/,$//')

	echo index_rem is ${index_rem}

	# Remove those that arent in bamlist cut complement removes them
	cat ${Sample}.beagle | cut --complement -f${index_rem} > ${Sample}.subset.beagle
	mv ${Sample}.subset.beagle  ${Sample}.beagle 

	echo $(expr $( cat ${Sample}.beagle | head -1 | tr "\t" "\n" | uniq | wc -l) - 3) samples remaining after subsetting

else 
	echo "no samples to remove"
fi

# Clean up
rm header
rm search.txt
rm to_remove

#--------------------------------------------------------------------------------
#-                           		Prune SNPs                                  -
#--------------------------------------------------------------------------------

if [ -n "$sitelist" ]; then
	echo "Subsetting to only those sites within sitelist: ${sitelist}"
    cat $(basename ${sitelist} .gz) | tr ":" "_" > sites_to_keep
	cat ${Sample}.beagle | head -n 1  > ${Sample}.subset.beagle
	cat ${Sample}.beagle | grep -Fwf sites_to_keep >> ${Sample}.subset.beagle
	mv ${Sample}.subset.beagle ${Sample}.beagle
else 
	echo "no sites to remove"
fi

#--------------------------------------------------------------------------------
#-                                    Run NGSDIST                               -
#--------------------------------------------------------------------------------
# Zip the beagle file again
pigz ${Sample}.beagle -p ${SLURM_CPUS_PER_TASK} --fast

# Prep variables and files
nsites=$(pigz -cd ${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | tail -n +2 | wc -l) 
nind=$(pigz -cd ${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | head -1 | sed $'s/\t/\\\n/g' | uniq | grep -v  "marker" | grep -v  "allele" | wc -l)
echo ${nsites} sites
echo ${nind} samples

# Pos file
pigz -cd ${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | awk '{print $1}' | sed 's/\(.*\)_/\1\t/' | tail -n +2 > ${Sample}.pos

# Labels file
pigz -cd ${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | head -1 | cut --complement -f 1,2,3 | tr "\t" "\n" | uniq > ${Sample}.labels

# Get total sites for just the analysed chromosomes
awk -v FS="\t" -v OFS="\t" '{print $1, $2}' ${ReferenceGenome}.fai > chr_lengths.txt
cat ${Sample}.pos | awk -v FS="\t" -v OFS="\t" '{print $1}' | sort | uniq > used_chrs.txt
totsites=$(cat chr_lengths.txt | grep -Fwf used_chrs.txt | cut -f2 | paste -sd+ | bc)

mkdir results

ngsDist --geno ${Sample}.beagle.gz -probs --pos ${Sample}.pos --labels ${Sample}.labels \
--n_ind ${nind} --n_sites ${nsites} --pairwise_del \
--evol_model 0 --n_boot_rep 100 --boot_block_size 1000 \
--out ${outname}.dist --n_threads ${SLURM_CPUS_PER_TASK}

# --tot_sites ${totsites}

# Split Main distance file into seperate files per bootstrap replicate
cat ${outname}.dist | awk -v RS= '{print > ("boot" NR ".dist")}'
ls | grep 'boot[0-9]' > distlist.txt
while read i; do
	oldcount=$(echo $i | grep -o '[0-9]\+')
	# If oldcount = 1, rename it to all (as its the full non-bootstrapped distances)
	if [ "$oldcount" -eq 1 ]; then
		newname=$(echo $i | sed -e "s/boot1/${outname}_all/g")
	else
        newcount=$(expr $oldcount - 1)
        newname=$(echo $i | sed -e "s/boot/${outname}_boot/g" | sed -e "s/boot${oldcount}/boot${newcount}/g")
	fi
    # Write out the file - replacing header with sample names
    paste <(echo -e "sample") <(cat ${Sample}.labels | tr '\n' '\t') > results/${newname}
    cat ${i} | tail -n +2 >> results/${newname}
done <distlist.txt
 
#--------------------------------------------------------------------------------
#-                                   Make tree                                  -
#--------------------------------------------------------------------------------
# Infer a tree for each of the matrices using FASTME
fastme -i ${outname}.dist -s -D 101 -o results/${outname}.nwk -I --nb_threads=${SLURM_CPUS_PER_TASK}

#split the allsites tree from the bootstraped ones
head -n 1 results/${outname}.nwk > results/${outname}.all.nwk
tail -n +2 results/${outname}.nwk | awk 'NF' > results/${outname}.boot.nwk

# Place supports on tree using using raxml
raxml-ng --support --tree results/${outname}.all.nwk --bs-trees results/${outname}.boot.nwk --prefix results/${outname} --threads ${SLURM_CPUS_PER_TASK}

# Copy files back to drive
cp results/* ${outdir}/${outname}

# Output useful job stats
/usr/local/bin/showJobStats.scr 
