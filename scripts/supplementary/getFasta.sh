#!/bin/bash
#SBATCH --job-name=angsdFASTA       
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
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

bamlist=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo bamlist=${bamlist}

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
while getopts ":R:I:O:q:m:r" options; do       
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
    I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;
    O)
      output=${OPTARG}
      echo output=${output}
    ;;
	q)
      basequal=${OPTARG}
      echo basequal=${basequal}
    ;;
    m)
      mapqual=${OPTARG}
      echo mapqual=${mapqual}
    ;;
	r)
	  echo "-r has been set, removing any sites that are variable in ancestral" >&2
	  remove_snps='true'
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

# Get sample name from vcf file
Sample=$(basename ${bamlist} | sed 's/bamlist_//g' | sed 's/manifest_//g' | sed 's/.txt//g')
echo Sample=${Sample}


#set outdir file name
if [[ $sitelist ]]; then
    outname=$(echo $Sample $(basename ${sitelist} | sed 's/.sites.*$//g' ) | sed 's/ /-/g' )
else 
    outname=$(echo $Sample denovo | sed 's/ /-/g' )
fi

# Make directories for outputs
mkdir -p ${output}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS or crams for input job
find "$(/usr/bin/ls -d ${indir})" \
| grep -F -f "${bamlist}" \
| grep -E "/($(tr '\n' '|' < "${bamlist}" | sed 's/|$//'))\.(bam|cram)$" \
| sort -u \
> "${Sample}_tmp.txt"

# Create a list of just filenames (no extension)
sed 's!.*/!!' "${Sample}_tmp.txt" \
  | grep -E -v '\.(bai|crai)$' \
  | grep -E '\.(bam|cram)$' \
  > "${Sample}_bams.txt"
  
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_bams.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

# Copy data files to temp and decompress
cp ${bamlist} .
cp ${ReferenceGenome} .
sleep 5
cp ${ReferenceGenome}.fai .

export ReferenceGenome


#Load Modules
module purge
#module load angsd/20250306-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load parallel/20240722-GCCcore-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

#--------------------------------------------------------------------------------
#-                          Copy and subset bam files                           -
#--------------------------------------------------------------------------------

if [[ $sitelist ]]; then
	# Subset bam files to just those reads overlapping the target sites
	echo subsetting bams to only reads overlapping target regions

	# Check if the file is compressed
	if [[ "${sitelist}" =~ \.gz$ ]]; then
		# If the file is gzipped, use zcat (or gzip -dc or gunzip -c) to decompress it on the fly
		zcat ${sitelist} > sitelist.txt
	else
		# If the file is not gzipped, process it directly
		cat ${sitelist} > sitelist.txt
	fi

	# Check if the file matches BED or Sites format
	if grep -qE '^[^[:space:]:#][^:\t#]*:[0-9]+(\r)?$' sitelist.txt; then
		echo "Sites file"
		awk -F':' -v OFS='\t' '
			/^[[:space:]]*#/ || /^[[:space:]]*$/ { next }          # skip comments/blank
			{ sub(/\r$/, "", $0) }                                  # strip CR if CRLF
			NF==2 && $2 ~ /^[0-9]+$/ { print $1, $2-1, $2 }         # 0-based half-open bed
		' sitelist.txt > sites.bed
	# Else, check for 3-column BED (and make awk fail if none)
	elif awk -F'\t' '
			BEGIN{ ok=0 }
			/^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
			NF==3 { ok=1 }
			END{ exit ok?0:1 }
		' sitelist.txt
	then
		echo "BED file sitelist"
		# If CRLF is possible, normalize while copying
		awk '{ sub(/\r$/, "", $0); print }' sitelist.txt > sites.bed

	else
		echo "Unknown format"
	fi
	echo Sites file contains $(bedtools makewindows -b sites.bed -w 1 | wc -l) sites
	
	# Run samtools view in parallel to subset and copy across, then index
	cat "${Sample}_tmp.txt" | parallel -j "${SLURM_CPUS_PER_TASK}" '
	  in={}
	  base=$(basename "$in")
	  stem=${base%.*}          # filename without last extension
	  ext=${base##*.}          # bam or cram

	  out="./${stem}.bam"
	  samtools view -b -T "'"${ReferenceGenome}"'" -L sites.bed "$in" > "$out"
	  samtools index "$out"
		
	  echo "subset ${stem} -> ${out}"
	'
else
	
	# Copy/convert files across, then index
	parallel -j "${SLURM_CPUS_PER_TASK}" --linebuffer --halt now,fail=1 --env ReferenceGenome '
	  in="{}"
	  base=$(basename -- "$in")
	  stem=${base%.*}          # filename without last extension
	  ext=${base##*.}          # bam or cram

	  out="./${stem}.bam"

	  samtools view -b -T "$ReferenceGenome" -- "$in" > "$out"
	  samtools index -- "$out"

	  echo "copied ${stem} -> ${out}"
	' :::: "${Sample}_tmp.txt"
fi

find "$PWD" -maxdepth 1 -type f -name '*.bam' -print > "${Sample}_bams.txt"

# Check BAMs
samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'

#--------------------------------------------------------------------------------
#-                                GetFASTA                                      -
#--------------------------------------------------------------------------------

# Keep only sites present in half of samples
nind=$(cat ${Sample}_bams.txt | wc -l)
nind=$(( nind / 2 ))

# launch angsd
if [[ "$remove_snps" = true ]]; then
	echo "Removing all sites that are variable"
	~/angsd/angsd/angsd -bam ${Sample}_bams.txt \
		-ref $(basename ${ReferenceGenome}) \
		-minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
		-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
		-doMaf 1 -doMajorMinor 1 \
		--ignore-RG 0 \
		-doCounts 1 -SNP_pval 0.01 -rmSNPs 1 -minMaf 0 \
		-GL 2 -minind ${nind} \
		-nThreads ${SLURM_CPUS_PER_TASK} \
		-doFasta 2 \
		-explode 1 \
		-out ${outname} 		
else  
	echo "Keeping all variable sites but selecting most common allele"
	~/angsd/angsd/angsd -bam ${Sample}_bams.txt \
		-ref $(basename ${ReferenceGenome}) \
		-minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
		-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
		-doMaf 1 -doMajorMinor 1 \
		--ignore-RG 0 \
		-doCounts 1 -SNP_pval 1.0 -minMaf 0 \
		-GL 2 -minind ${nind} \
		-nThreads ${SLURM_CPUS_PER_TASK} \
		-doFasta 2 \
		-explode 1 \
		-out ${outname} 

fi

# Unzip fasta and index
pigz -d -p ${SLURM_CPUS_PER_TASK} ${outname}.fa.gz
samtools faidx ${outname}.fa

module purge
module load seqtk/1.4-GCC-13.3.0

# Count number of N bases in output
total_bases=$(seqtk comp ${outname}.fa | awk '{sum += $2} END {print sum}')
n_bases=$(seqtk comp ${outname}.fa | awk '{sum += $9} END {print sum}')
proportion_n=$(echo "scale=4; $n_bases / $total_bases" | bc)
echo "$n_bases of $total_bases bases in ${outname}.fa are N ($proportion_n %)"

# Copy files back to drive
cp -r ${outname}.fa ${output}/.
cp -r ${outname}.fa.fai ${output}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
