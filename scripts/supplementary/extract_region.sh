#!/bin/bash
#SBATCH --job-name=ext_region         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=72:00:00
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

# This script searches the genome for the coordinates of all sequences in a reference fasta, or a bed coordinates file
# then extracts those regions from BAM files and calls variants

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=extraction_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi


#--------------------------------------------------------------------------------
#-                          Parse Input Arguments                               -
#--------------------------------------------------------------------------------

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R Reference Genome ] [ -B Bed file ] [ -F Fasta file ] [ -O Input directory ] [ -O Output directory ] [ -l minimum length ] [ -d minimum depth ] [ -p padding on each side of region ] [ -d minimum depth ] [ -t ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1

while getopts ":R:B:F:I:O:d:l:p" options; do       
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
    B)
      bedfile=${OPTARG}
      if [[ $bedfile ]]; then
	  # Test if exists
        if [ ! -f "$bedfile" ] ; then  
          echo "Error: -B ${bedfile} doesnt exist"
          exit_abnormal
          exit 1
        fi
      fi
	  echo bedfile=${bedfile}	  
      ;;
    F)                             
      ref_fasta=${ref_fasta}
      if [[ $ref_fasta ]]; then
          # Test if exists
          if [ ! -f "$ref_fasta" ] ; then  
            echo "Error: -F ${ref_fasta} doesnt exist"
            exit_abnormal
            exit 1
          fi
      fi
	  echo ref_fasta=${ref_fasta}	 
    ;;
    I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;        
    O)
      outdir=${OPTARG}
      echo outdir=${outdir}
      ;;
    d)
      mindepth=${OPTARG}
      if [[ $mindepth ]]; then
        echo mindepth=${mindepth}
      else
        mindepth=1
        echo using default mindepth of 1 read
      fi
    ;; 
    l)
      minlength=${OPTARG}
      if [[ $minlength ]]; then
        echo minlength=${minlength}
      else
        minlength=20
        echo using default minlength of 20bp
      fi
    ;; 
    p)
      padding=${OPTARG}
      if [[ $padding ]]; then
        echo padding=${padding}
      else
        padding=0
        echo using default padding of 0bp
      fi
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
input=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
Sample=$(basename ${input} .txt | sed -r 's/bamlist_//')

# Make directories for outputs
mkdir ${outdir}

# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Process manifest file
# This must be a tab delimited file with format:
# Col1 = sample id
# Col2 (optional) = population/species (Considered seperately in variant calling)
# Col3 = output name for each sample
ncol=$(awk -v 'FS=\t' '{print NF}' ${input} | sort -nu | tail -n 1)
if [[ "$ncol" -eq 2 ]]; then
    echo 'No group column provided in input, considering all samples as single population for variant calling'
    cat ${input} | sed 's/^.*\///g'| sed 's/.bam//g' > samples.txt
    cat ${input} | awk '{ print $1, $2}' | sed 's/^.*\///g'| sed 's/.bam//g' > new_names.txt
elif [[ "$ncol" -eq 3 ]]; then
    echo 'Group column provided in input, considering each group as separate population for variant calling'

    cat ${input} | awk '{ print $1}' | sed 's/^.*\///g'| sed 's/.bam//g' > samples.txt
    cat ${input} | awk -v 'FS=\t' -v 'OFS=\t' '{ print $1, $2}' | sed 's/ /_/g' | sort | uniq > groups.txt
    cat ${input} | awk '{ print $1, $3}' | sed 's/^.*\///g'| sed 's/.bam//g' > new_names.txt
else 
    echo 'input file must only have 2 or 3 columns'
    exit
fi

## Create list of BAM files to process
find $(/usr/bin/ls -d ${indir}) | grep -F -f samples.txt | sort | uniq  > ${Sample}_tmp.txt
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_bams.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"


# Copy all data files to temp directory and decompress
cp $(dirname ${ReferenceGenome})/* .
cp ${ReferenceGenome}.fai .
xargs -a ${Sample}_tmp.txt cp -ft .
pigz -p ${SLURM_CPUS_PER_TASK} -d ./*.gz

# Load modules
module purge
module load BLAT/3.5-GCC-8.2.0-2.31.1
module load bedops/2.4.35
module load SAMtools/1.17-GCC-11.2.0
module load BEDTools/2.30.0-GCC-11.2.0
module load MAFFT/7.453-gompi-2020a-with-extensions
module load BCFtools/1.17-GCC-11.2.0

#--------------------------------------------------------------------------------
#-                         	  Set up regions file                         		-
#--------------------------------------------------------------------------------

# Check if reference fasta was provided
rm to_extract.bed
touch to_extract.bed
if [[ $ref_fasta ]]; then
    # conduct blat search for the regions in fasta
    blat $(basename ${ReferenceGenome}) ${ref_fasta} alignments.psl

    # convert results to BED format
    psl2bed < alignments.psl >> to_extract.bed
fi

# Check if BED file of coordinates was provided
if [[ $bedfile ]]; then
    cat $bedfile >> to_extract.bed
fi


# add extra padding bases to each side and replace any negative numbers with 0
awk -v padding=${padding} 'BEGIN { FS="\t"; OFS=FS } {$2-=padding;$3+=padding}1' to_extract.bed | sed 's/-[0-9]\+/1/g' > tmpbed
mv tmpbed to_extract.bed

# make a regions file from the bed file for SAMtools
awk 'BEGIN { FS="\t" } { print $1":"$2"-"$3}' to_extract.bed > to_extract.regions

# Create list of bams to process
ls | grep .bam$ > bamlist.txt

#--------------------------------------------------------------------------------
#-                        Call variants within regions                    		-
#--------------------------------------------------------------------------------

# Check if groups are defined in the manifest file
if [ ! -f groups.txt ] ; then 
    # If true - Call variants, grouping samples by population/species
    bcftools mpileup -f ${ReferenceGenome} -b bamlist.txt -R to_extract.bed -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP,FORMAT/QS | \
    bcftools call -mv -G groups.txt  | \
    bcftools norm -f ${ReferenceGenome} -Ou |\
    bcftools filter --IndelGap 5 -Ou |\
    bcftools sort - -Oz -o ${Sample}.vcf.gz

else 
    # If false - Call variants considering all samples as the same population
    bcftools mpileup -f ${ReferenceGenome} -b bamlist.txt -R to_extract.bed -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP,FORMAT/QS | \
    bcftools call -mv | \
    bcftools norm -f ${ReferenceGenome} -Ou |\
    bcftools filter --IndelGap 5 -Ou |\
    bcftools sort - -Oz -o ${Sample}.vcf.gz
fi

# Index resulting VCF
tabix ${Sample}.vcf.gz 

#--------------------------------------------------------------------------------
#-                    	  Create consensus fasta files                          -
#--------------------------------------------------------------------------------
# Call consensus FASTA file for each sample from VCF file
mkdir fastas
# Loop through each sample id
for i in $(bcftools query -l ${Sample}.vcf.gz );	do 
    #Masking approach modified from https://github.com/samtools/bcftools/issues/1386
    echo ${i}
    
    # Subset joint vcf file to just that individual and retain variant sites only 
    bcftools view ${Sample}.vcf.gz -m2 -M2 -s ${i} --min-ac 1 -Ob -o tmp.bcf
    tabix tmp.bcf
    
    # Subset bam to target regions to count coverage
    samtools view -L to_extract.bed -b ${i}.bam > subset.bam
    
    # Find low coverage sites (those below mindepth) and output bedgraph format - these sites will be masked with N bases
    bedtools genomecov -bga -ibam subset.bam | awk -v m=${mindepth} '$4 < m' > low_coverage_sites.bed
    cat low_coverage_sites.bed > mask.bed
    
    #Commented out code below will include variants, even if they are at low coverage
    #bcftools query -f'%CHROM\t%POS0\t%END\t%QUAL\n' tmp.bcf > variants.bed
    #bedtools subtract -a low_coverage_sites.bed -b variants.bed > mask.bed
    
    # Loop through regions
    while read l; do
      echo ${l}
      # Subset regions file to just that marker
      echo ${l} > tmp.regions
     
      # Call consensus for regions in the bed file - Masking low coverage sites
      samtools faidx ${ReferenceGenome} -r tmp.regions  | bcftools consensus tmp.bcf -s ${i} -p ${i}_ --mark-del '-' -m mask.bed -H I  > fastas/${i}_${l}.fa
    done <to_extract.regions

done

#--------------------------------------------------------------------------------
#-                    	    Filter and Align fastas                             -
#--------------------------------------------------------------------------------

# Merge all seperate FASTA files by loci, applying lengthfilter, then align with MAFFT
# Loop through regions
while read l; do
  outname=$(echo ${l} "dp"${mindepth} | sed 's/ /_/g' | sed 's/:/_/g')
  rm ${outname}.fa

  # Extract the target region from the reference genome to align against
  echo ${l} | tr ':' '\t' | tr '-' '\t' > tmp.bed
  bedtools getfasta -fi ${ReferenceGenome} -bed tmp.bed > ${outname}.fa
  
  # Get a list of fasta files for that region
  find fastas | grep ${l} > tmpfastas
  
  # Filter records to remove those with < lengthfilter bases
  while read f; do
    total_bases=$(grep -v ">" ${f} | wc | awk '{print $3-$1}')
	n_bases=$(grep -v ">" ${f} | tr -cd 'N|n' | wc -c)
	diff_bases=$(expr $total_bases - $n_bases)	
	if [[ "$diff_bases" -gt "$minlength" ]]; then
        # Rename sample using column 3 of the input file
        current_name=$(cat ${f} | grep '>' | sed 's/_.*$//g' | sed 's/>//g')
        new_name=$(grep "$current_name" new_names.txt | awk '{ print $2}' FS=' ')
        # Passes lengthfilter, output record to multi-sample fasta
        cat ${f} | sed "s/$current_name/$new_name/g" >> ${outname}.fa
	else
        # Doesnt past lengthfilter, exclude record
        echo less than ${minlength} bases in ${f}
	fi    
  done <tmpfastas
  
  # check if any of the samples passed lengthfilter (i.e. FASTA file has records)
  if grep -q '>' "${outname}.fa" ; then
    # Align records with MAFFT
    mafft ${outname}.fa > ${outname}_aligned.fa
    cp ${outname}_aligned.fa ${outdir}/${outname}.fa
  else
    # No records passed lengthfilter
    echo no data in  ${outname}.fa
  fi
done <to_extract.regions

# Output useful job stats
/usr/local/bin/showJobStats.scr 