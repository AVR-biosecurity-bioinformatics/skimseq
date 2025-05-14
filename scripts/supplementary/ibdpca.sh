#!/bin/bash
#SBATCH --job-name=ibdpca         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --time=240:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                        Install sofware from github                           -
#--------------------------------------------------------------------------------

# Install specific version of angsd where IBS was working
#cd ~
#module purge
#module load HTSlib/1.21-GCC-13.3.0
#git clone https://github.com/angsd/angsd.git;
#cd angsd
#git reset --hard 68b0838
#make

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

Index=ibdpca_job_index.txt

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
ibstype=1

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
while getopts ":R:S:I:O:q:m:ct" options; do       
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
	c)
	  echo "-c has been set using consensus instead of single read sampling" >&2
	  ibstype=2
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

# Get sample name from vcf file
Sample=$(basename ${bamlist} | sed 's/bamlist_//g' | sed 's/.txt//g')
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


## Create list of BAMS for input job
find "$(/usr/bin/ls -d ${indir})" | grep -F -f "${bamlist}" | grep -E "/($(tr '\n' '|' < ${bamlist} | sed 's/|$//')).bam$" | sort | uniq > "${Sample}_tmp.txt"
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_bams.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"


# Copy data files to temp and decompress
cp ${bamlist} .
cp ${ReferenceGenome} .
sleep 5
cp ${ReferenceGenome}.fai .

#Load Modules
module purge
#module load angsd/20250306-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load parallel/20240722-GCCcore-13.3.0

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
	if grep -P '^[A-Za-z0-9\.]+:\d+$' "sitelist.txt" > /dev/null; then
		echo "Sites file"
		# Processing Sites format assuming each line is like CM028320.1:2010752
		awk -F':' '{print $1, $2-1, $2}' sitelist.txt > sites.bed
	# Check if the file matches BED format
	elif awk -F'\t' 'NF == 3' sitelist.txt > /dev/null; then
		echo "BED file sitelist"
		# Assuming sitelist.txt is already in BED format with 3 fields
		cp sitelist.txt sites.bed
	else
		echo "Unknown format"
	fi
	echo Sites file contains $(bedtools makewindows -b sites.bed -w 1 | wc -l) sites
	
	# Run samtools view in parallel to subset and copy across, then index
	cat ${Sample}_tmp.txt | parallel -j ${SLURM_CPUS_PER_TASK} "samtools view -b -L sites.bed {} > ./{/.}.bam && samtools index ./{/.}.bam && echo subset {/.}"

else
	# Copy files across, then index
	cat ${Sample}_tmp.txt | parallel -j ${SLURM_CPUS_PER_TASK} "cp {} . && samtools index ./{/.}.bam && echo copied {/.}"
fi
#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------

# if -t is set, all BAM files will be subset to the first chr and only 10% of data
if [[ "$test_run" = true ]]; then
	fractionOfReads=0.1
	chr1=$( samtools idxstats $(ls | grep '.bam$' | grep -v '.bai' | head -n 1) | head -1 | awk -F '\t' '{print $1}')
	echo subsampling BAMs to only ${chr1} and ${fractionOfReads} of reads for testing 
		
for i in $(cat ${Sample}_bams.txt);do
    echo $i
        samtools view -@ ${SLURM_CPUS_PER_TASK} -b1 $i ${chr1} -s $fractionOfReads > output.bam
        mv output.bam $i
        samtools index $i -@ ${SLURM_CPUS_PER_TASK}
    done
fi

# Check BAMs
samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'

#--------------------------------------------------------------------------------
#-                                      Run ANGSD                               -
#--------------------------------------------------------------------------------

if [[ $sitelist ]]; then
    echo "only analysing sites from ${sitelist}"
	
	# Expand bedfile all sites and transform to angsd sitelist
	bedtools makewindows -b sites.bed -w 1 | awk -v OFS="\t" -v FS="\t" '{print $1, $3}' > sites.txt

	# Index sites file
	~/angsd/angsd/angsd sites index sites.txt

	# Using -sites to subset to target sites
	~/angsd/angsd/angsd -bam ${Sample}_bams.txt -sites sites.txt \
		-ref $(basename ${ReferenceGenome}) \
		-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
		-minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
		-GL 2 -doMajorMinor 1 \
		-doIBS ${ibstype} -doCounts 1 -doCov 1 -makeMatrix 1 \
		-nThreads ${SLURM_CPUS_PER_TASK} \
		-out ${outname} 
else 
    echo 'Sitelist not provided, estimating sites denovo'
	# Using -sites to subset to target sites
	
	~/angsd/angsd/angsd -bam ${Sample}_bams.txt \
		-ref $(basename ${ReferenceGenome}) \
		-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
		-minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
		-GL 2 -doMajorMinor 1 \
		-doIBS ${ibstype} -doCounts 1 -doCov 1 -makeMatrix 1 \
		-nThreads ${SLURM_CPUS_PER_TASK} \
		-out ${outname} 
fi

# Removed these parameters
# -doMaf 1 -minMaf 0.01

# Count non-missing data per individual
zcat ${outname}.ibs.gz | awk 'NR>1 {for (i=5; i<=NF; i++) count[i] += ($i != "-1")} END {for (i=5; i<=NF; i++) print count[i]}' > tmp
paste <(cat ${Sample}_bams.txt) <(cat tmp) > ${outname}.indcounts.txt

# Count non-missing data per site
zcat ${sitelist} > sites.txt
zcat ${outname}.ibs.gz | tail -n+2 > genotypes.txt

awk '
BEGIN { FS=OFS="\t" }
NR==FNR { sites[$1] = NR; order[NR] = $1; next } # Read second_file into array with order tracking
NR>1 {
    site = $1 ":" $2
    count = 0
    for (i=5; i<=NF; i++) count += ($i != "-1")
    output[site] = $1 "\t" $2 "\t" $3 "\t" $4 "\t" count
}
END {
    for (i=1; i<=length(order); i++) {
        site = order[i]
        if (site in output) {
            print output[site]
        } else {
            split(site, parts, ":")
            print parts[1], parts[2], "NA", "NA", "NA"
        }
    }
}' sites.txt genotypes.txt > ${outname}.sitecounts.txt

# Copy files back to drive
cp ${outname}*.covMat ${output}/${outname}/.
cp ${outname}*.ibs.gz ${output}/${outname}/.
cp ${outname}*.ibsMat ${output}/${outname}/.
cp ${Sample}_bams.txt ${output}/${outname}/.
cp ${outname}.indcounts.txt ${output}/${outname}/.
cp ${outname}.sitecounts.txt ${output}/${outname}/.

rm -rf *
# Output useful job stats
/usr/local/bin/showJobStats.scr 
