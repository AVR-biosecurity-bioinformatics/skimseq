#!/bin/bash
#SBATCH --job-name=trim_align_merge         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20gb
#SBATCH --time=72:00:00
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

# This script trims reads, aligns paired end fastq reads to a genome, then merges bams

#--------------------------------------------------------------------------------
#-                                  Parse SBATCH                                -
#--------------------------------------------------------------------------------

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

# Get sample info from submission
FullSampleName=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
Sample=$(basename $FullSampleName | sed 's/_S[0-9].*$//')  # | cut -d'_' -f1
SequencePath=$( find $(/usr/bin/ls -d ${SLURM_SUBMIT_DIR}/fastq) -maxdepth 3 -name '*.fastq.gz' -type f | grep [\/_-]${Sample}_ | sort | uniq )

echo FullSampleName=${FullSampleName}
echo Sample=${Sample}
echo SequencePath=${SequencePath}

[[ ! -z "${Sample}" ]] && echo $(echo ${SequencePath} | wc -w | awk '{print $1/2}') file/s to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
ReferenceGenome=""
MitoGenome=""
known_variants=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R ReferenceGenome ] [ -M MitoGenome ] [ -V known_variants ] [ -O otutput directory ] [ -c ] [ -i ] [ -u ] [ -q ] [ -t ]" 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":R:M:V:O:icuqt" options; do       
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
	  align_genome="true"
	  echo ReferenceGenome=${ReferenceGenome}	  
      ;;
    M)
      MitoGenome=${OPTARG}
	  # Test if exists
	  if [ ! -f "$MitoGenome" ] ; then  
        echo "Error: -M ${MitoGenome} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  align_mito="true"
	  echo MitoGenome=${MitoGenome}
      ;;
	V)
	  known_variants=${OPTARG}
	  # Test if not empty
	  if [[ $known_variants ]]; then
		  # test if exists
		  if [ ! -f "$known_variants" ] ; then  
			echo "Error: -V ${known_variants} doesnt exist"
			exit_abnormal
			exit 1
		  fi
		echo "Base qualities will be recalibrated using known variants" >&2  
		recalibrate="true"
		echo known_variants=${known_variants}
	  fi
      ;;
    O)
      outdir=${OPTARG}
      echo outdir=${outdir}
    ;;
	i)
	  echo "Indels will be realigned using GATK IndelRealigner" >&2
	  realign_indels='true'
      ;;
	c)
	  echo "Variants will be called per-sample using GATK HaplotypeCaller" >&2
	  call_variants='true'
    ;; 
	u)
	  echo "unaligned reads will be extracted" >&2
	  extract_unaligned='true'
    ;;
	q)
	  echo "FASTQC will be run" >&2
	  run_fastqc='true'
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

# Check if mitogenome and reference genomes are indexed - if not exit with warning

#--------------------------------------------------------------------------------
#-                                Prepare files                                 -
#--------------------------------------------------------------------------------
# Make directories for outputs
mkdir -p ${outdir}
mkdir ${outdir}/qc
mkdir ${outdir}/mito
mkdir ${outdir}/offtarget
mkdir ${outdir}/bams
mkdir ${outdir}/gvcf

# Go to a temp directory to conduct analyses
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd $tmp_dir
pwd

# Copy data files to temp and decompress
echo ${SequencePath} | tr " " "\n" > ${Sample}_tmp.txt

cat ${Sample}_tmp.txt
xargs -a ${Sample}_tmp.txt cp -ft .
cp ${ReferenceGenome}* .
cp ${ReferenceGenome}.fai .
cp $(find $(dirname $ReferenceGenome) -name *dict -type f) .
cp ${MitoGenome}* .
pigz -p${SLURM_CPUS_PER_TASK} -d ./*.fastq.gz

ls

# Create read tracking file
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep '_R1_' | sed -r 's/.fastq.gz//' | sort | uniq  > ${Sample}_R1.txt
cat ${Sample}_R1.txt


#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------
# if -t is set, all fastqs will be subsampled to make the run fastQC
if [[ "$test_run" = true ]]; then
	echo "TEST RUN - SUBSAMPLING ALL FASTQS"
	module purge
	module load BBMap/39.17-GCC-13.3.0
	
	# Loop across fastq files
	for R1 in $(cat ${Sample}_R1.txt ) ;do 
	R2=$(echo ${R1} | sed -r 's/_R1_/_R2_/g')                                                                                               
	echo ${R1}
	echo ${R2}
	# Run bbmap reformat
	reformat.sh in1=${R1}.fastq in2=${R2}.fastq out1=${R1}_subsampled.fastq out2=${R2}_subsampled.fastq reads=1000 sampleseed=666 usejni=t
	mv ${R1}_subsampled.fastq ${R1}.fastq
	mv ${R2}_subsampled.fastq ${R2}.fastq
	done
fi
#--------------------------------------------------------------------------------
#-                                  FastQC	                                    -
#--------------------------------------------------------------------------------
# This step uses the fastqc tool to analyse sequence qualities

if [[ "$run_fastqc" = true ]]; then
	
  # Load modules
  module purge
  module load FastQC/0.12.1-Java-11
  
  # Run FastQC on all fastq files
  mkdir fastqc_pretrim
  fastqc *.fastq --outdir fastqc_pretrim
  
  # Copy fastQC results back to the qc folder in the origin location
  cp -r fastqc_pretrim ${outdir}/qc/.

else
	echo Skipping pre trimming FASTQC
fi

#--------------------------------------------------------------------------------
#-                                 Fastp trimming                               -
#--------------------------------------------------------------------------------
# This step uses the fastp tool to quality filter reads, and remove adapters and polyG tails
# Load modules
module purge
module load fastp/0.23.4-GCC-13.3.0

mkdir fastp
# Loop across fastq files
for R1 in $(cat ${Sample}_R1.txt ) ;do 
	R2=$(echo ${R1} | sed -r 's/_R1_/_R2_/g')                                                                                               
	echo ${R1}
	echo ${R2}

	# Run fastp
    fastp -i ${R1}.fastq -I ${R2}.fastq \
    -q 20 \
    --length_required 15 \
    --n_base_limit 5 \
    --trim_poly_g \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 20 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --correction \
    --overlap_len_require 30 \
    --overlap_diff_limit 5 \
    --overlap_diff_percent_limit 20 \
    --thread ${SLURM_CPUS_PER_TASK} \
    -o ${R1}.trimmed.fastq -O ${R2}.trimmed.fastq \
    -h fastp/${R1}.fastp.html -j fastp/${R1}.fastp.json -R ${Sample} 
    
    # Remove original fastqs from tempdir
    rm ${R1}.fastq 
    rm ${R2}.fastq 
done

# Copy fastp results back to origin location
cp -r fastp ${outdir}/qc/.

#--------------------------------------------------------------------------------
#-                           FastQC post trim	                                  -
#--------------------------------------------------------------------------------
# This step uses the fastqc tool to analyse sequence qualities

if [[ "$run_fastqc" = true ]]; then
	
  # Load modules
  module purge
  module load FastQC/0.12.1-Java-11
  
  # Run FastQC on all fastq files
  mkdir fastqc_posttrim
  fastqc *.trimmed.fastq --outdir fastqc_posttrim
  
  # Copy fastQC results back to the qc folder in the origin location
  cp -r fastqc_posttrim ${outdir}/qc/.

else
	echo Skipping post trimming FASTQC
fi

#--------------------------------------------------------------------------------
#-                             Extract mitochondria                             -
#--------------------------------------------------------------------------------
# This step aligns all reads to the mitochondrial genome and returns a consensus fasta file

if [[ "$align_mito" = true ]]; then
	echo Aligning reads to Mitochondrial genome $(basename "${MitoGenome}")

	## Load modules
	module purge
	module load BWA/0.7.18-GCCcore-13.3.0
	module load SAMtools/1.21-GCC-13.3.0
	module load BCFtools/1.21-GCC-13.3.0

	# Loop across fastq files
	for R1 in $(cat ${Sample}_R1.txt ) ;do 
		R2=$(echo ${R1} | sed -r 's/_R1_/_R2_/g')                                                                                               
		echo ${R1}
		echo ${R2}
		
		# Setup read group headers, these are necessary for merging of replicates
		RG_ID=$(echo ${R1} | awk -F _ '{print $1 "." $4}')
		RG_LB=$(echo ${R1} | awk -F _ '{print $2}')
		RG_PI=$(grep peak fastp/${R1}.fastp.json | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]')
		#RG_UID=$(echo ${R1} | awk -F _ '{print $1 "." $4 "_" $3}')
		
		# Align to Mitochondrial Genome using the bwa-mem algorithm
		bwa mem -t ${SLURM_CPUS_PER_TASK} \
		-R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${Sample}\tPI:${RG_PI}") \
		"$(basename "${MitoGenome}")" \
		${R1}.trimmed.fastq \
		${R2}.trimmed.fastq \
		| samtools view -bS - > ${Sample}_${R1}.tempmito.bam 
	done

	# Merge any replicate BAM files, and remove PCR duplicates
	find . -name "*.tempmito.bam" > to_merge.txt
	nbam=$(cat to_merge.txt | wc -l)

	# Check if to_merge.txt is larger than one
	if [[ $nbam -gt 1 ]]; then
		echo Merging $nbam BAM files
		samtools merge -b to_merge.txt -O bam - \
		| samtools sort -n -O BAM \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - ${Sample}.mito.bam
	else
		echo only one BAM file, skipping merge
		samtools sort -n ${Sample}_${R1}.tempmito.bam  -O bam \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - ${Sample}.mito.bam
	fi

	# Index resulting bam file
	samtools index -@ ${SLURM_CPUS_PER_TASK} ${Sample}.mito.bam

	# Call consensus sequence and output as fasta file
	bcftools mpileup -Ou -f ${MitoGenome} ${Sample}.mito.bam -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP | \
    bcftools call -Ou -m --ploidy 1 | \
    bcftools norm -f ${MitoGenome} -Ou |\
    bcftools filter --IndelGap 5 -Oz -o ${Sample}.vcf.gz
	tabix ${Sample}.vcf.gz 
    
	bcftools consensus -f ${MitoGenome} ${Sample}.vcf.gz --absent N --missing N  > ${Sample}.mito.fa
    
	# Copy just the fasta files back to origin
	cp ${Sample}.mito.fa ${outdir}/mito/.
    rm ${Sample}.vcf.gz 
    rm ${Sample}.mito.bam

	echo Mitochondria successfully aligned to $(basename "${MitoGenome}")

else
	echo no mitogenome file provided, skipping mitochondrial extraction
	# TODO add denovo mito assembly
fi

#--------------------------------------------------------------------------------
#-                              Align to genome                                 -
#--------------------------------------------------------------------------------
# This step aligns all reads to the reference genome

if [[ "$align_genome" = true ]]; then
	echo aligning reads to reference genome $(basename "${ReferenceGenome}")

	# Load modules
	module purge
	module load BWA/0.7.18-GCCcore-13.3.0
	module load SAMtools/1.21-GCC-13.3.0

	for R1 in $(cat ${Sample}_R1.txt ) ;do 
		R2=$(echo ${R1} | sed -r 's/_R1_/_R2_/g')                                                                                               
		echo ${R1}
		echo ${R2}
		
		# Setup read group headers, these are necessary for merging of replicates
		RG_ID=$(echo ${R1} | awk -F _ '{print $1 "." $4}')
		RG_LB=$(echo ${R1} | awk -F _ '{print $2}')
		RG_PI=$(grep peak fastp/${R1}.fastp.json | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]')
		#RG_UID=$(echo ${R1} | awk -F _ '{print $1 "." $4 "_" $3}')

		# Align reads to reference genome
		bwa mem -t ${SLURM_CPUS_PER_TASK} \
		-R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${Sample}\tPI:${RG_PI}") \
		"$(basename "${ReferenceGenome}")" \
		${R1}.trimmed.fastq \
		${R2}.trimmed.fastq \
		| samtools view -bS - > ${Sample}_${R1}.temp.bam 
        
        # Remove trimmed fastqs from tempdir
        rm ${R1}.trimmed.fastq
        rm ${R2}.trimmed.fastq
	done

	# Merge any replicate BAM files, and remove PCR duplicates
	find . -name "*.temp.bam" > to_merge.txt
	nbam=$(cat to_merge.txt | wc -l)

	# Check if to_merge.txt is larger than one
	if [[ $nbam -gt 1 ]]; then
		echo Merging $nbam BAM files
		samtools merge -b to_merge.txt -O bam - \
		| samtools sort -n -O BAM \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - ${Sample}.merged.bam
        
        # Remove temporary bam files
        cat to_merge.txt | xargs rm
	else
		echo only one BAM file, skipping merge
		samtools sort -n ${Sample}_${R1}.temp.bam  -O bam \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - ${Sample}.merged.bam
        
        # Remove temporary bam files
        rm ${Sample}_${R1}.temp.bam        
	fi		
    
	# Index resulting bam file
	samtools index -@ ${SLURM_CPUS_PER_TASK} ${Sample}.merged.bam ${Sample}.merged.bam.bai

	# Make sure all BAMs are properly formatted
	samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'
	echo Reads successfully aligned to $(basename "${ReferenceGenome}")

else
	echo no reference genome file provided, exiting script
	exit 1
fi

#--------------------------------------------------------------------------------
#-                        	  Extract unaligned reads		                    -
#--------------------------------------------------------------------------------
# This step extracts all unaligned reads into fastq files for manual exploration
if [[ "$extract_unaligned" = true ]]; then
	echo extracting unaligned reads

	# Extract unaligned reads
	samtools bam2fq -@ ${SLURM_CPUS_PER_TASK} -1 ${Sample}.unmapped.R1.fastq.gz -2 ${Sample}.unmapped.R2.fastq.gz -0 /dev/null -s /dev/null -f12 ${Sample}.merged.bam

	# Move to output directory
	mv ${Sample}.unmapped.R1.fastq.gz ${outdir}/offtarget/.
	mv ${Sample}.unmapped.R2.fastq.gz ${outdir}/offtarget/.

	# Reindex bam file
	samtools index -@ ${SLURM_CPUS_PER_TASK} ${Sample}.merged.bam 
	echo unaligned reads extracted
else
	echo Skipping unalgined read extraction
fi

#--------------------------------------------------------------------------------
#-                        Realign INDELs using GATK realigner                   -
#--------------------------------------------------------------------------------
# This step realigns reads around any detected INDEL mutations to ensure they arent alignment artefacts
# TODO: Validate this is required 

if [[ "$realign_indels" = true ]]; then
	echo realigning reads around indels

	# Load earlier gatk which has indel realigner
	module purge
	module load GATK/3.8-0-Java-1.8

	# Create a list of target sites for realignment
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
		-R $(basename ${ReferenceGenome}) \
		-I ${Sample}.merged.bam \
		-o ${Sample}.merged.intervals

	# Realign around target sites
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
		-R $(basename ${ReferenceGenome}) \
		-I ${Sample}.merged.bam \
		-targetIntervals  ${Sample}.merged.intervals \
		-o ${Sample}.realigned.bam

	echo Indel realignment complete
else
	echo Skipping unalgined read extraction
	mv ${Sample}.merged.bam ${Sample}.realigned.bam
fi

# Remove merged bam
rm ${Sample}.merged.bam

#--------------------------------------------------------------------------------
#-                 		   Base Quality Score Recalibration                     -
#--------------------------------------------------------------------------------
# This step recalibrates base quality scores and is critical for merging across multiple sequencing runs
# This step requires a list of known variants, if these are not provided the step will be skipped 

# Check if a list of known variants was provided
if [[ "$recalibrate" = true ]]; then
	echo Using known variants from ${known_variants} to recalibrate base qualities
	module purge
	module load GATK/4.6.1.0-GCCcore-13.3.0-Java-21
	module load R/4.4.2-gfbf-2024a
	
	# Transform sitelist into bed file
	pigz -cd ${known_variants} -p ${SLURM_CPUS_PER_TASK} | awk -v FS=':' -v OFS='\t' '{print $1,$2-1,$2}' > known_variants.bed
	
	# Index bed file
	gatk IndexFeatureFile \
		-I known_variants.bed

	# Create recalibration table
	gatk BaseRecalibrator \
		-I ${Sample}.realigned.bam \
		-R $(basename ${ReferenceGenome}) \
		-known-sites known_variants.bed \
		-O ${Sample}.recal_data.table1 

	# Run BQSR to recalibrate quality scores
	gatk ApplyBQSR \
		-R $(basename ${ReferenceGenome}) \
		-I ${Sample}.realigned.bam \
		-bqsr ${Sample}.recal_data.table1  \
		-O ${Sample}.recal.bam 


	# Analyse variation in the adjusted BAM
	gatk BaseRecalibrator \
		-I ${Sample}.recal.bam \
		-R $(basename ${ReferenceGenome}) \
		-known-sites known_variants.bed \
		-O ${Sample}.recal_data.table2 

	# make before and after plot
	gatk AnalyzeCovariates \
		-before ${Sample}.recal_data.table1 \
		-after ${Sample}.recal_data.table2 \
		-plots ${Sample}.recalplots.pdf

	mv ${Sample}.recal.bam ${Sample}.bam 
    
    mkdir recal
    cp ${Sample}.recal_data.table* recal/.
	cp ${Sample}.recalplots.pdf recal/.
    cp -r recal ${outdir}/qc/.

else
	echo no known_variants file provided skipping base quality recalibration
	mv ${Sample}.realigned.bam ${Sample}.bam 
fi

# Remove realigned bam
rm ${Sample}.realigned.bam

#--------------------------------------------------------------------------------
#-                 			Generate BAM stats 			   			            -
#--------------------------------------------------------------------------------
# This final step ensures the BAM is indexed properly and outputs coverage statistics
module purge
module load SAMtools/1.21-GCC-13.3.0

# Index final BAM
samtools index -@ ${SLURM_CPUS_PER_TASK} ${Sample}.bam ${Sample}.bam.bai

# Output sample coverage statistics
mkdir stats
samtools coverage ${Sample}.bam > stats/${Sample}.stats.txt
cp -r stats ${outdir}/qc/.

# Copy BAM files back to drive
cp ${Sample}.bam ${outdir}/bams/.
cp ${Sample}.bam.bai ${outdir}/bams/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
