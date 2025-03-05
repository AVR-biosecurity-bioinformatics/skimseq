#!/bin/bash
#SBATCH --job-name=relate         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=6
#SBATCH --mem=60GB 
#SBATCH --time=640:00:00
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

# This script uses NGSRelate to estimate relatedness between samples

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=relate_job_index.txt

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
while getopts ":B:S:M:O:t" options; do       
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
    M)                             
      maf=${OPTARG}
	  # Test if exists
	  if [ ! -f "$maf" ] ; then  
        echo "Error: -M ${maf} doesnt exist"
        exit_abnormal
        exit 1
      fi
      freqs="true"
	  echo maf=${maf}	  
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
Sample=$(basename ${beagle} .beagle.gz)
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

# Copy data files to temp and decompress
cp ${beagle} .
cp ${sitelist} .
cp ${bamlist} .

ls 
pigz -p ${SLURM_CPUS_PER_TASK} -d *.gz

#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------
# if -t is set, all fastqs will be subsampled to make the run fastQC
if [[ "$test_run" = true ]]; then
    # For quick testing - subsample input to just first chromosome
    echo subsampling to first chr/contig only
    chr1=$(cat ${Sample}.beagle | head -n 2 | tail -n 1 | awk -F '\t' '{print $1}' | awk -F '_' '{print $1}')
    cat ${Sample}.beagle | head -n 1  > ${Sample}.test.beagle
    cat ${Sample}.beagle | grep ${chr1} >> ${Sample}.test.beagle
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
#-                                 Run ngsRelate                                -
#--------------------------------------------------------------------------------
# Prep variables and files
nsites=$(cat ${Sample}.beagle | tail -n +2 | wc -l) 
nind=$(cat ${Sample}.beagle | head -1 | sed $'s/\t/\\\n/g' | uniq | grep -v  "marker" | grep -v  "allele" | wc -l)

echo ${nsites} sites
echo ${nind} samples

#Load Modules
module purge
module load NgsRelate/2.0-GCC-11.2.0

pigz ${Sample}.beagle -p ${SLURM_CPUS_PER_TASK} --fast

if [[ "$freqs" = true ]]; then
    echo "Using allele frequencies file"
    cat $(basename ${sitelist} .gz) | tr ":" "\t" > sites_to_keep
    cp ${maf} .
    pigz -d $(basename ${maf})
    awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1 OFS $2]; next}
        ($1 OFS $2) in a {
            print $7
        }
        ' sites_to_keep $(basename ${maf} .gz) >> freq
     
    # Run NGSRelate
    ngsRelate \
    -G ${Sample}.beagle.gz \
    -f freq \
    -n ${nind} \
    -L ${nsites} \
    -p ${SLURM_CPUS_PER_TASK} \
    -O ${outname}.relate \
    -l 0.05
else
    echo "Estimating allele frequencies from beagle file"
    # Run NGSRelate
    ngsRelate \
    -G ${Sample}.beagle.gz \
    -n ${nind} \
    -L ${nsites} \
    -p ${SLURM_CPUS_PER_TASK} \
    -O ${outname}.relate \
    -l 0.05
fi
# Rename output columns to the actual sample name using join (requires sorting on the join columns)

# Create sample index file for renaming relate output
zcat ${Sample}.beagle.gz | head -1 | sed $'s/\t/\\\n/g' | uniq | grep -v  "marker" | grep -v  "allele" | awk -v OFS='\t' '{ print NR-1, $1 }' | sort -k 1 > sample_index.txt

cat ${outname}.relate | sort -k 2 >  ${outname}.relate.sorted 

#Update the second field first using join
join sample_index.txt ${outname}.relate.sorted -1 1 -2 2 -t $'\t' | cut -f2- | sort -k2 > tmp.txt

#Update the second field using join
cat ${outname}.relate | head -1 > ${outname}.relate.renamed 
join sample_index.txt tmp.txt -1 1 -2 2 -t $'\t' | cut -f2- >> ${outname}.relate.renamed 
mv ${outname}.relate.renamed  ${outname}.relate

# Copy files back to drive
cp ${outname}.relate ${outdir}/${outname}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
