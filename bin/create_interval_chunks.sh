#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = ref_genome
# $4 = counts_per_chunk
# $5 = split_overweight
# $6 = min_interval_gap
# $7 = included_contigs
# $8 = include_zero
# $9 = counts_file

# This script uses bed files with or without a counts column to assign bed intervals to seperate chunks for parallel processing
# For Haplotypecaller, input is a single bed with per-base counts in column 4
# For GenotypeGVCFS, input is multiple bed files with no counts
# For Mpileup, input is multiple bed fils with per-base counts in column 4

TARGET_COUNTS_PER_CHUNK=$(awk -v x="${4}" 'BEGIN {printf("%d\n",x)}')
OUTDIR=$(pwd)
SPLIT_OVERWEIGHT="${5}"
GAP_BP="${6}"
ALL_BASES="${8}"
TMPDIR=$(mktemp -d)
CPUS="${1}"

# Create contig bed file including just those in the included intervals file
awk 'NR==FNR{
        if($0 !~ /^#/ && $1!="") c[$1]=1
        next
     }
     ($1 in c){ print $1, 0, $2 }' OFS=$'\t' ${7} ${3}.fai > contigs.bed
     

# Exit early if contigs.bed is empty
if [[ ! -s contigs.bed ]]; then
  : > "_empty.bed"          # create empty output file
  exit 0
fi

btmp="${TMPDIR}/tmp.bed"

# Extract just included contigs from zipped bed hten lexagraphically sort all counts files
parallel -j "${CPUS}" --line-buffer '
  f={}
  out="${f%.bed.gz}.sorted.bed"
  tabix "$f" -R contigs.bed | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$out"
' :::: counts_files.list

# Merge pre-sorted files with bedops, then sum across GAP_BP, sort back into genome (can be non lexagraphical) order.
bedops --everything *.sorted.bed \
| bedtools merge -i - -d ${GAP_BP} -c 4 -o sum \
| bedtools sort -i - -g ${3}.fai > "$btmp"

# If ALL_BASES=false: keep only intervals that have any counts 
# If ALL_BASES=true: build whole-contig intervals, map counts back, drop zeros
if [[ "$ALL_BASES" == "false" ]]; then
  mv "$btmp" intervals_with_counts.bed
else
  bedtools map -sorted -g ${3}.fai -a contigs.bed -b "$btmp" -c 4 -o sum -null 0 \
  | awk -v OFS='\t' '$4 != 0' \
  > intervals_with_counts.bed
fi

# Exit early if intervals_with_counts.bed is empty
if [[ ! -s intervals_with_counts.bed ]]; then
  : > "_empty.bed"          # create empty output file
  exit 0
fi

# Split intervals that individually exceed the target counts.
# Assumes counts are roughly uniform across the interval length.
if [[ "$SPLIT_OVERWEIGHT" == "true" ]]; then
  awk -v target="$TARGET_COUNTS_PER_CHUNK" '
  BEGIN{OFS="\t"}
  {
    chrom=$1; start=$2; end=$3; w=$4
    len=end-start

    # guard against weird zero/negative lengths
    if (len <= 0) next

    # if not overweight, keep as-is
    if (w <= target) {
      print chrom, start, end, w
      next
    }

    # number of pieces needed (ceil)
    n = int((w + target - 1) / target)

    # don’t create more pieces than bases
    if (n > len) n = len

    base = int(len / n)
    rem  = len - base * n

    substart = start
    for (i=1; i<=n; i++) {
      sz = base + (i <= rem ? 1 : 0)
      subend = substart + sz

      # proportional weight by length to preserve density
      subw = w * sz / len

      # round to integer
      subw = int(subw + 0.5)
      print chrom, substart, subend, subw
      substart = subend
    }
  }
  ' intervals_with_counts.bed > intervals_split.bed
else
  cat intervals_with_counts.bed > intervals_split.bed
fi

# Calculate total counts and number of intervals
TOTAL_COUNTS=$( awk '{s+=$4} END{print s+0}' intervals_split.bed )
N_INTERVALS=$(cat intervals_split.bed | wc -l)

# Decide number of file splits (chunks) to keep counts balanced.
K=$(( (TOTAL_COUNTS + TARGET_COUNTS_PER_CHUNK - 1) / TARGET_COUNTS_PER_CHUNK )) 
if [ "$K" -lt 1 ]; then K=1; fi
if [ "$K" -gt "$N_INTERVALS" ]; then K="$N_INTERVALS"; fi

# Assign intervals to chunks, keep it contiguous so they are in sorted order
awk -v outdir="$OUTDIR" -v tot="$TOTAL_COUNTS" -v n="$N_INTERVALS" -v K="$K" '
function abs(x){ return x<0 ? -x : x }

BEGIN{
  chunk=1
  remK=K
  remW=tot
  sum=0
  ideal = remW / remK
  fname = sprintf("%s/chunk_%d.bed", outdir, chunk)
}

{
  chrom=$1; start=$2; end=$3; w=$4+0
  lines_after = n - FNR

  if (sum>0 && remK>1) {
    stop_now_better = (abs(sum-ideal) <= abs((sum+w)-ideal))
    enough_left     = (lines_after >= (remK-2))

    if (stop_now_better && enough_left) {
      chunk++
      remK--
      sum=0
      ideal = remW / remK
      fname = sprintf("%s/chunk_%d.bed", outdir, chunk)
    }
  }

  print chrom"\t"start"\t"end"\t"w > fname
  sum  += w
  remW -= w
}
' intervals_split.bed
   
# Rename each output file
for i in *chunk_*.bed;do
  # Pad output chunk names
  n=$(basename "$i" | sed -E 's/^chunk_([0-9]+)\.bed/\1/')
  pad=$(printf "%05d" "$n")

  # compute hash of this chunk’s contents
  hash=$(md5sum "$i" | awk '{print $1}')

  # Output just column 1:4
  out="_${pad}${hash}.bed"
  cut -f1-4 "$i" > "$out"
  
  # report size
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$i")
  counts=$(awk '{s+=$4} END{print s+0}' "$i")
  echo "${out}: ${bases} genomic bases, ${counts} summed counts"

  rm $i
done
