#!/bin/bash
#SBATCH --job-name=abbababa       
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=400GB
#SBATCH --time=72:00:00
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

Index=abbababa_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

manifest=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo manifest=${manifest}

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
while getopts ":R:A:S:I:O:r" options; do       
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
	  # Test if exists
	  if [ ! -f "$AncestralGenome" ] ; then  
        echo "Error: -A ${AncestralGenome} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo AncestralGenome=${AncestralGenome}	  
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
    r)
	  echo "-r has been set, using regions file (slow for lots of sites)" >&2
	  use_regions='true'
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
# Create output name
outname=$( basename $manifest | sed 's/manifest_//g' |sed 's/\.txt//g'   )

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS for input job
find $(/usr/bin/ls -d ${indir}) | grep -F -f ${manifest} | sort | uniq  > ${outname}_tmp.txt
cat ${outname}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${outname}_bams.txt
[[ ! -z "${outname}" ]] && echo $(wc -l ${outname}_tmp.txt | awk '{ print$1 }') BAM files to process for ${outname} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

# Copy data files to temp and decompress
cp ${ReferenceGenome} .
sleep 20
cp ${ReferenceGenome}.fai .
cp ${AncestralGenome} .
sleep 20
cp ${AncestralGenome}.fai .
xargs -a ${outname}_tmp.txt cp -ft .

#Load Modules
module purge
module load GSL/2.7-GCC-11.2.0
module load SAMtools/1.15-GCC-11.2.0
module load HTSlib/1.15-GCC-11.2.0
module load BCFtools/1.15-GCC-11.2.0
module load R/4.2.0-foss-2021b
#module load angsd/0.933-GCC-9.3.0

# Check BAMs
samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'

#--------------------------------------------------------------------------------
#-                                ABBA BABA                                     -
#--------------------------------------------------------------------------------


if [[ $sitelist ]]; then
    echo "only analysing sites from ${sitelist}"

    # Check if sitelist is tab delimited (augmented) or colon delimited (simple)
    tmp=$(zcat ${sitelist} | head -1)
    if echo "$tmp" | grep -q ':' ; then
        echo 'Sitelist is semicolon delimited - assuming its a simple sites file'
        # Setup sitelist if provided and run angsd
        pigz -cd ${sitelist} -p ${SLURM_CPUS_PER_TASK} | tr ':' '\t' | sort -V -k1,1 -k2,2n > sites.txt
    elif echo "$tmp" | grep -q '[[:space:]]' ; then
        echo 'Sitelist is tab delimited - assuming its an augmented sites file'
        pigz -cd ${sitelist} -p ${SLURM_CPUS_PER_TASK} | sort -V -k1,1 -k2,2n > sites.txt
    fi

    # Index sites file
    /home/ap0y/angsd/angsd/angsd sites index sites.txt

    # REGIONS IS TAKING UP LONGEST
    if [[ "$use_regions" = true ]]; then
        # Just keep the chromosome and positions column for next step
        cat sites.txt | awk -v FS='\t' -v OFS='\t' '{{ print $1,$2}}' > tmpsites

        # Split regions file every 1000 SNPS - per chromosome
        echo -n "" > rf_starts.txt
        cut -f1 sites.txt | sort | uniq > chrs.txt
        while read c; do
            cat tmpsites | grep $c | head -1 > tmp.txt
            cat tmpsites | grep $c | sed -n '0~1000p' >> tmp.txt
            cat tmpsites | grep $c | tail -1 >> tmp.txt
            cat tmp.txt | uniq >> rf_starts.txt
        done <chrs.txt 
        nlines=$(cat rf_starts.txt | wc -l)

        echo -n "" > rf.txt
        for i in $(seq 1 $(cat rf_starts.txt | wc -l)); do
         nextline=$(($i + 1))
         start_chr=$(cat rf_starts.txt | sed -n ${i}p | cut -f1) 
         end_chr=$(cat rf_starts.txt | sed -n ${nextline}p | cut -f1) 
         if [ "$start_chr" = "$end_chr" ]; then
             start_pos=$(cat rf_starts.txt | sed -n ${i}p | cut -f2) 
             end_pos=$(($(cat rf_starts.txt | sed -n ${nextline}p | cut -f2) -1))
             echo "$start_chr":"$start_pos"-"$end_pos" >> rf.txt
         else
            # Check if its the final line, or just a chromosome break
            if [ "$i" -eq "$nlines" ]; then
                # last line - dont write anything
                echo 'last line in file'
            else
                echo 'moving to next chromosome'
            fi
         fi
        done

        cat rf.txt
        cat rf.txt | sort > rf2.txt 

        # Calculate ABBA BABA using regions and sites file
        /home/ap0y/angsd/angsd/angsd -bam ${outname}_bams.txt -sites sites.txt -rf rf2.txt \
        -ref $(basename ${ReferenceGenome}) -anc $(basename ${AncestralGenome}) \
        -minMapQ 20 -baq 2 -C 50 -minQ 15 \
        -remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
        -doCounts 1 \
        -doAbbababa 1 \
        -nThreads ${SLURM_CPUS_PER_TASK} \
        -out ${outname} 
    else
        # Calculate ABBA BABA - restricting to the sitelist but not using regions
        /home/ap0y/angsd/angsd/angsd -bam ${outname}_bams.txt -sites sites.txt -rf rf2.txt \
        -ref $(basename ${ReferenceGenome}) -anc $(basename ${AncestralGenome}) \
        -minMapQ 20 -baq 2 -C 50 -minQ 15 \
        -remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
        -doCounts 1 \
        -doAbbababa 1 \
        -nThreads ${SLURM_CPUS_PER_TASK} \
        -out ${outname} 
    fi
     
     
else
    echo 'Sitelist not provided, estimating sites denovo'
    /home/ap0y/angsd/angsd/angsd -bam ${outname}_bams.txt \
    -ref $(basename ${ReferenceGenome}) -anc $(basename ${AncestralGenome}) \
    -minMapQ 20 -baq 2 -C 50 -minQ 15 \
    -remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
    -doCounts 1 \
    -doAbbababa 1 \
    -nThreads ${SLURM_CPUS_PER_TASK} \
    -out ${outname} 

fi

#estimate Z score
Rscript /home/ap0y/angsd/angsd/R/jackKnife.R file=${outname}.abbababa indNames=${outname}_bams.txt outfile=${outname}_abbababa

# Copy files back to drive
cp -r ${outname}_abbababa.txt ${output}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
