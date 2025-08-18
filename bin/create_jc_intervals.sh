#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = jc_interval_scaling_factor
# $4 = include_bed     
# $5 = exclude_bed
# $6 = Reference_genome
# $7 = jc_interval_min_n

# NOTE: Intervals for joint calling are made directly from the chromosomes
# any masked bases should not have been included in single sample calling

# ---- GATK Resources ----
# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# ---- Calculate nchunks ----
# Number of samples
n_samples=$(ls *.g.vcf.gz | wc -l)

# Count the number of records in the largest gvcf (by filesize)
n_variants=$(bcftools view -H $(ls -S *.g.vcf.gz | head -n1) | wc -l)

available_mem_gb=48          # Mem available for the joint calling steps (in GB)
safety_factor=0.8            # fraction of memory to actually use
scaling_factor=$(awk -v x="${3}" 'BEGIN {printf("%.10f\n", x)}')  # GB per variant per sample (2.5e-7) - ROUGH ESTIMATE

# Total memory required if full genome in one interval
mem_total=$(echo "$n_variants * $n_samples * $scaling_factor" | bc -l)

# Max memory per chunk considering safety factor
max_mem_per_chunk=$(echo "$available_mem_gb * $safety_factor" | bc -l)

# Recommended number of chunks (round up)
n_chunks=$(echo "($mem_total / $max_mem_per_chunk + 0.9999)/1" | bc)
if (( n_chunks < ${7} )); then n_chunks=${7}; fi

# ---- Create chromosome list ----
# Create an included intervals file, only contig names that are in this file will be retained
bedtools subtract -a ${4} -b ${5} > included_intervals.bed
awk '{contigs[$1]=1} END {for(c in contigs) print c}' included_intervals.bed > included_contigs.txt

# Calculate number of bases that should theoretically be contained in each interval
exp_bases_per_group=$(awk -v n_chunks="$n_chunks" '{sum += $2} END {print sum/n_chunks}' ${6}.fai)

# Subset chromosomes longer than expected bases per group, and filter to just those contig names in included intervals
awk -v minlen="$exp_bases_per_group" '{ if($2 >= minlen) print $1 "\t0\t" $2 }' ${6}.fai > long.bed
awk 'NR==FNR {keep[$1]; next} ($1 in keep)' included_contigs.txt long.bed > long.filtered.bed

# Subset chromosomes shorter than expected bases per group, and filter to just those contig names in included intervals
awk -v minlen="$exp_bases_per_group" '{ if($2 < minlen) print $1 "\t0\t" $2 }'  ${6}.fai > short.bed
awk 'NR==FNR {keep[$1]; next} ($1 in keep)' included_contigs.txt short.bed > short.filtered.bed

mv long.filtered.bed long.bed
mv short.filtered.bed short.bed

# Total bases in long contigs - these contigs will be split into smaller intervals
long_bases=$(awk '{sum += $3 - $2} END {print sum}' long.bed)

# Total bases in short contigs - these contigs will be grouped together
short_bases=$(awk '{sum += $3 - $2} END {print sum}' short.bed)

# Total bases in genome
total_bases=$((long_bases + short_bases))

# Divide the total number of splits between the long and short sets
long_splits=$(( n_chunks * long_bases / total_bases ))
short_splits=$(( n_chunks * short_bases / total_bases ))

# Ensure at least 1 interval per set, as the short contigs go through different processing to long ones
if (( long_splits < 1 )); then long_splits=1; fi
if (( short_splits < 1 )); then short_splits=1; fi

# Divide the long chromosomes into shorter ones, DONT MIX CONTIGS
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SplitIntervals \
   -R ${6} \
   -L long.bed \
   --dont-mix-contigs true \
   --scatter-count ${long_splits} \
   --interval-merging-rule OVERLAPPING_ONLY \
   --extension -long.interval_list \
   -O $(pwd) \
   --subdivision-mode INTERVAL_SUBDIVISION

# Group the short chromosomes together, MIX CONTIGS
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SplitIntervals \
   -R ${6} \
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

