#!/bin/bash
#SBATCH --job-name=extract         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB 
#SBATCH --time=24:00:00
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

Index=extract_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

# Get sample info from submission
bedfile=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
vcf=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -B Bedfile ] [ -O Outdir ] [ -t ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":B:O:t" options; do       
  # use silent error checking;
  case "${options}" in
    V)                             
      vcf=${OPTARG}
	  # Test if exists
	  if [ ! -f "$vcf" ] ; then  
        echo "Error: -B ${vcf} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo vcf=${vcf}	  
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

# Make directories for outputs
mkdir ${outdir}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Copy data files to temp and decompress
cp ${vcf} tmp.vcf.gz
cp ${vcf}.tbi tmp.vcf.gz.tbi
cp ${bedfile} targets.bed

# Load modules
module purge
module load bedops/2.4.35
module load SAMtools/1.15-GCC-11.2.0
module load BCFtools/1.15-GCC-11.2.0
module load MAFFT/7.453-gompi-2020a-with-extensions

#--------------------------------------------------------------------------------
#-                            	  Extract regions                         		-
#--------------------------------------------------------------------------------

mkdir fasta
while read r; do
  # Create regions file 
  echo ${r} > tmp.bed
  awk 'BEGIN { FS=" " } { print $1":"$2"-"$3}' tmp.bed > tmp.regions
  
  # Create output filename
  outname=$(awk 'BEGIN { FS=" " } { print $4}' tmp.bed)
  rm ${outname}.fa
  touch ${outname}.fa
  for i in $(bcftools query -l tmp.vcf.gz );	do 
      # Filter vcf by sample - keeping only pass tags
      bcftools view tmp.vcf.gz -s ${i} -r $(cat tmp.regions) -Oz -o tmp2.vcf.gz
      tabix -f -p vcf tmp2.vcf.gz

      # Check if any did 
      n_var=$(zcat tmp2.vcf.gz | grep -v "^#" | wc -l)
      echo "$n_var variants remain in locus ${outname} for sample ${i}"
      if [ "$n_var" -gt 0 ]; then
        # Output IUPAC characters for each allele
        #samtools faidx /group/referencedata/mspd-db/genomes/insect/daktulosphaera_vitifoliae/V3.1/Dv_genome_V3.1.fa -r tmp.regions | \
        #bcftools consensus tmp2.vcf.gz -s ${i} -H I > tmp.fa
        #cat tmp.fa | sed -r "s/>/>${i}_${outname}_/" >> ${outname}.fa
      
        # Output first phased allele and rename headers
        samtools faidx /group/referencedata/mspd-db/genomes/insect/daktulosphaera_vitifoliae/V3.1/Dv_genome_V3.1.fa -r tmp.regions | \
        bcftools consensus tmp2.vcf.gz -s ${i} -H 1 > tmp.fa
        cat tmp.fa | sed -r "s/>/>${i}_${outname}_1_/" > allele1.fa
        
        
        # Output second phased allele and rename headers
        samtools faidx /group/referencedata/mspd-db/genomes/insect/daktulosphaera_vitifoliae/V3.1/Dv_genome_V3.1.fa -r tmp.regions  | \
        bcftools consensus tmp2.vcf.gz -s ${i} -H 2 > tmp.fa
        cat tmp.fa | sed -r "s/>/>${i}_${outname}_2_/" > allele2.fa
        
        cat allele1.fa  >> ${outname}.fa
        cat allele2.fa  >> ${outname}.fa
      fi
  done
  cp ${outname}.fa fasta/.
  rm ${outname}.fa
  # Align outputs
  #mafft --maxiterate 1000 --localpair ${outname}.fa > ${outname}_aligned.fa
  #mv ${outname}_aligned.fa fasta/.
  #rm ${outname}.fa
done <targets.bed

cp -r fasta/* ${outdir}

# Output useful job stats
/usr/local/bin/showJobStats.scr 