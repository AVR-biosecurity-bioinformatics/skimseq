#!/bin/bash
#SBATCH --job-name=GATK_joint         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB 
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

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=gatk_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

interval=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo interval=${interval}

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
while getopts ":R:M:I:O:t" options; do       
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
    M)
      manifest=${OPTARG}
      echo manifest=${manifest}
      ;;
	I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;
    O)
      output=${OPTARG}
      echo output=${output}
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

# Exit if no ref genome supplied
if [[ "$ReferenceGenome" = "" ]]; then
	exit_abnormal
fi

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------

Sample=$(basename ${output})
echo ${Sample}

# Create output name
n_int=$(echo ${interval} | tr -d -c ',' | awk '{ print length; }')
if [ ! -z "$n_int" ]; then
    first_chr=$(echo ${interval} | cut -d "," -f1 )
    last_chr=$(echo ${interval} | rev | cut -d "," -f1 | rev )
    outname=$(echo ${Sample}_${first_chr}-${last_chr})
else
    outname=$(echo ${Sample}_${interval})
fi
echo ${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd
find $(/usr/bin/ls -d ${indir}) | grep -F -f ${manifest} | sort | uniq > ${Sample}_tmp.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_tmp.txt | awk '{ print$1 }') VCF files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

dictfile=$(find $(dirname $ReferenceGenome) -maxdepth 1 -name *dict -type f ) 

# Copy data files to temp and decompress
cp ${ReferenceGenome} .
cp ${ReferenceGenome}.fai .
cp ${dictfile} .
xargs -a ${Sample}_tmp.txt cp -ft .
ls -l

##Create cohort mapping file 
#needs to look like this:
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gz

find $(pwd) -name "*.g.vcf.gz"  | grep -v  ".tbi" | sort | uniq  > ${Sample}_gvcfs.txt
cat ${Sample}_gvcfs.txt | sed 's=.*/==' | sed 's/.g.vcf.gz//g' > ${Sample}_names.txt
paste ${Sample}_names.txt ${Sample}_gvcfs.txt > cohort.sample_map

#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------
# if -t is set, all fastqs will be subsampled to make the run fastQC
if [[ "$test_run" = true ]]; then
    echo "TEST RUN - Subsetting to first chr for testing"

    # Load BCFtools for filtering (assuming module loading is required)
    module purge
    module load BCFtools/1.9-intel-2019a
    
    # Find name of first interval
    first_chr=$(echo "${interval}" | cut -d "," -f1)

    # Loop through VCF files listed in ${Sample}_gvcfs.txt
    while read vcf_file; do
        echo "${vcf_file}"
        
        # Subset to the first chromosome for quick testing
        temp_vcf="${vcf_file/.g.vcf.gz/_subset.g.vcf.gz}"
        bcftools view "${vcf_file}" --regions "${first_chr}" -Oz -l 1 -o "${temp_vcf}"
        mv "${temp_vcf}" "${vcf_file}"
        
        # Reindex the subsetted VCF
        tabix -p vcf "${vcf_file}"
    done < "${Sample}_gvcfs.txt"
fi

#--------------------------------------------------------------------------------
#-                                      Run GATK                               -
#--------------------------------------------------------------------------------
# Load modules 
module purge
module load shifter/22.02.1

# Create bed file with target intervals
echo ${interval} | tr ',' '\n' > intervals.txt
cat intervals.txt | tr ':' '\t' | tr '-' '\t' > intervals.bed

# convert bed to interval list for gatk
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk BedToIntervalList -I intervals.bed  -O $(basename ${ReferenceGenome}).interval_list -SD $(basename ${dictfile})

# Make genomic DB from gvcf files
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk --java-options "-Xmx45G" GenomicsDBImport \
	  --genomicsdb-workspace-path ${outname}.db \
	  --sample-name-map cohort.sample_map \
	  --batch-size 50 \
	  --tmp-dir ${tmp_dir} \
	  --reader-threads ${SLURM_CPUS_PER_TASK} \
	  --genomicsdb-shared-posixfs-optimizations true \
	  -L $(basename ${ReferenceGenome}).interval_list
    
#Joing genotype gvcfs
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk --java-options "-Xmx45G" GenotypeGVCFs \
	  -R $(basename ${ReferenceGenome}) \
	  -V gendb://${outname}.db \
	  -O ${outname}.vcf.gz \
	  --genomicsdb-shared-posixfs-optimizations true \
	  --tmp-dir ${tmp_dir} 
 
# Copy files back to drive
cp ${outname}.vcf.gz*  ${output}

echo completed GATK run

# Output useful job stats
/usr/local/bin/showJobStats.scr 
