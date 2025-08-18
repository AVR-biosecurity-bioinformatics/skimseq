#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = interval_n
# $4 = Reference_genome

# NOTE: Intervals for joint calling are made directly from the chromosomes, any masked bases should not have been included in single sample calling
# This is to properly handle genomicsDBimport 

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# Convert any scientific notation to integers
interval_n=$(awk -v x="${3}" 'BEGIN {printf("%d\n",x)}')

# Calculate number of bases that should theoretically be contained in each interval
exp_bases_per_group=$(awk -v interval_n="$interval_n" '{sum += $2} END {print sum/interval_n}' ${4}.fai)

# Subset chromosomes longer than expected bases per group - these will be split into smaller intervals
awk -v minlen="$exp_bases_per_group" '{ if($2 >= minlen) print $1 "\t0\t" $2 }'  ${4}.fai > long.bed

# Subset chromosomes shorter than expected bases per group - these will be grouped together
awk -v minlen="$exp_bases_per_group" '{ if($2 < minlen) print $1 "\t0\t" $2 }'  ${4}.fai > short.bed

# Total bases in long contigs
long_bases=$(awk '{sum += $3 - $2} END {print sum}' long.bed)

# Total bases in short contigs
short_bases=$(awk '{sum += $3 - $2} END {print sum}' short.bed)

# Total bases in genome
total_bases=$((long_bases + short_bases))

# Number of splits for long and short sets
long_splits=$(( interval_n * long_bases / total_bases ))
short_splits=$(( interval_n * short_bases / total_bases ))

# Ensure at least 1 interval per set
if (( long_splits < 1 )); then long_splits=1; fi
if (( short_splits < 1 )); then short_splits=1; fi

# Divide the long chromosomes into shorter ones, DONT MIX CONTIGS
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SplitIntervals \
   -R ${4} \
   -L long.bed \
   --dont-mix-contigs true \
   --scatter-count ${long_splits} \
   --interval-merging-rule OVERLAPPING_ONLY \
   --extension -long.interval_list \
   -O $(pwd) \
   --subdivision-mode INTERVAL_SUBDIVISION

# Group the short chromosomes together, MIX CONTIGS
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SplitIntervals \
   -R ${4} \
   -L short.bed \
   --scatter-count ${short_splits} \
   --interval-merging-rule OVERLAPPING_ONLY \
   --extension -short.interval_list \
   -O $(pwd) \
   --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW

# Rename and convert each split output into a bed, and add padding
for i in *.interval_list;do
  # Hashed output name
  HASH=$( md5sum "$i" | awk '{print $1}' ) 

  # Convert resulting interval list to bed format
  java "-Xmx${java_mem}G" -jar $EBROOTPICARD/picard.jar IntervalListToBed \
  	--INPUT $i \
  	--OUTPUT tmp.bed
	
	#adding extra character at start to ensure that other temp beds dont accidentally get passed to next process
    cut -f1-4 "tmp.bed" > "_${HASH}.bed"
	
  # remove intermediate files
  rm $i tmp.bed
done

# Remove temporary bed files
rm -f long.bed short.bed

