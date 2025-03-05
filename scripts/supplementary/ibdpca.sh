#!/bin/bash
#SBATCH --job-name=ibdpca         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --time=240:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
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

Index=ibdpca_job_index.txt

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
while getopts ":R:S:O:ct" options; do       
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
    O)
      output=${OPTARG}
      echo output=${output}
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
# Read input params
bamlist=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})

# Get sample name from vcf file
Sample=$(basename ${bamlist} | sed 's/bamlist_//g' | sed 's/.txt//g')
echo Sample=${Sample}

#set output file name
outname=$(echo $Sample $(echo  $(basename ${sitelist} | sed 's/.sites.*$//g' )) | sed 's/ /-/g' )

# Make directories for outputs
mkdir -p ${output}/${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS for input job
find $(/usr/bin/ls -d ${SLURM_SUBMIT_DIR}/bams/*) | grep -F -f ${bamlist} | grep "\.bam$" | sort | uniq  > ${Sample}_tmp.txt
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_bams.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Copy reference files to temp and decompress
cp ${ReferenceGenome} .
sleep 5
cp ${ReferenceGenome}.fai .

#Load Modules
module purge
module load SAMtools/1.15-GCC-11.2.0
module load HTSlib/1.15-GCC-11.2.0
module load BCFtools/1.15-GCC-11.2.0
module load parallel/20210722-GCCcore-11.2.0

#--------------------------------------------------------------------------------
#-                              Subset bam files                                -
#--------------------------------------------------------------------------------

# Subset bam files to just those reads overlapping the target sites
echo subsetting bams to only reads overlapping target regions

# Create bedfile 
zcat ${sitelist} | awk -F: '{OFS="\t"; print $1, $2-1, $2}' > sites.bed

# Run samtools view in parallel to subset and copy across
cat ${Sample}_tmp.txt | parallel -j ${SLURM_CPUS_PER_TASK} "samtools view -b -f 1 -L sites.bed {} > ./{/.}.bam && samtools index ./{/.}.bam && echo subset {/.}"

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

# Setup sitelist
zcat ${sitelist} | tr ':' '\t' | sort -V -k1,1 -k2,2n > sites.txt

##Old regions file code - doesnt seem to improve performance
#cut -f1 sites.txt |sort | uniq > chrs.txt
#
# SPlit regions file every 1000 SNPS - per chromosome
#echo -n "" > rf_starts.txt
#
#while read c; do
#    cat sites.txt | grep $c | head -1 > tmp.txt
#    cat sites.txt | grep $c | sed -n '0~1000p' >> tmp.txt
#    cat sites.txt | grep $c | tail -1 >> tmp.txt
#    cat tmp.txt | uniq >> rf_starts.txt
#done <chrs.txt 
#nlines=$(cat rf_starts.txt | wc -l)
#
#echo -n "" > rf.txt
#for i in $(seq 1 $(cat rf_starts.txt | wc -l)); do
# nextline=$(($i + 1))
# start_chr=$(cat rf_starts.txt | sed -n ${i}p | cut -f1) 
# end_chr=$(cat rf_starts.txt | sed -n ${nextline}p | cut -f1) 
# if [ "$start_chr" = "$end_chr" ]; then
#     start_pos=$(cat rf_starts.txt | sed -n ${i}p | cut -f2) 
#     end_pos=$(($(cat rf_starts.txt | sed -n ${nextline}p | cut -f2) -1))
#     echo "$start_chr":"$start_pos"-"$end_pos" >> rf.txt
# else
#    # Check if its the final line, or just a chromosome break
#    if [ "$i" -eq "$nlines" ]; then
#        # last line - dont write anything
#        echo 'last line in file'
#    else
#        echo 'moving to next chromosome'
#    fi
# fi
#done
#
#cat rf.txt
#cat rf.txt | sort > rf2.txt 


# Index sites file
/home/ap0y/angsd/angsd/angsd sites index sites.txt

# launch angsd (local github version) 
# Using -sites to subset to target sites
/home/ap0y/angsd/angsd/angsd -bam ${Sample}_bams.txt -sites sites.txt \
-ref $(basename ${ReferenceGenome}) \
-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
-minMapQ 20 -baq 2 -C 50 -minQ 15 \
-GL 2 -doMajorMinor 1 \
-doIBS ${ibstype} -doCounts 1 -doCov 1 -makeMatrix 1 \
-nThreads ${SLURM_CPUS_PER_TASK} \
-out ${outname} 

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
