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
SPLIT_OVERWEIGHT=${5}
GAP_BP="${6}"


# Flag to include zero counts (keeps entire contigs, used for short contig chunking for genotypegvcfs) or exclude
if [[ ${8} == "false" ]]; then
    BG="-bg"
else 
    BG="-bga"
fi

# Create contig list in reference order with lengths, including just those in the included intervals file
awk 'NR==FNR{c[$1]=1; next} c[$1]' ${7} ${3}.fai \
| cut -f1,2 > contigs.tsv

# Function for determining how many bed files (samples) overlap each position of the reference genome, and optionally counting bases
process_contig() {
  chr="$1"
  len="$2"

  out="per_contig/${chr}.bed"
  tmp_bg="per_contig/${chr}.bg.tmp"
  all_bg="per_contig/${chr}.all_intervals.bg"
  all_bed="per_contig/${chr}.all_intervals.bed"

  genome_line=$(printf "%s\t%s\n" "$chr" "$len")

  # Build bedGraph of intervals from the counts files that overlap positions of the genome
  # BG toggles whether intervals with no bases should be inclucded
  for f in *counts.bed.gz; do
    tabix "$f" "$chr" 2>/dev/null || true
  done \
  | bedtools genomecov $BG -i - -g <(printf "%s" "$genome_line") \
  > "$tmp_bg"

  # No records on this contig -> emit empty and stop
  if [ ! -s "$tmp_bg" ]; then
    : > "$out"
    rm -f "$tmp_bg"
    return 0
  fi

  # Merge adjacent/nearby segments and sum the depth column
  bedtools merge -i "$tmp_bg" -d "$GAP_BP" -c 4 -o sum > "$all_bg"

  # Map the original count values (col 4) from ALL counts files onto the intervals
  cut -f1-3 "$all_bg" > "$all_bed"

  for f in *counts.bed.gz; do
    tabix "$f" "$chr" 2>/dev/null || true
  done \
  | bedtools sort -i - \
  | bedtools map -a "$all_bed" -b - -c 4 -o sum -null 0 \
  > "$out"

  rm -f "$tmp_bg" "$all_bg" "$all_bed"
}

mkdir -p per_contig
export -f process_contig
export GAP_BP
export MODE
export BG

parallel --colsep '\t' --jobs "${1}" process_contig {1} {2} :::: contigs.tsv

# Merge contig counts
: > combined_counts.bed
while IFS=$'\t' read -r chr len; do
  cat "per_contig/${chr}.bed" >> combined_counts.bed
done < contigs.tsv

# TODO: Could add filter to remove chunks with not many samples called before merging intervals

# Exit early if combined_counts.bed is empty
if [[ ! -s combined_counts.bed ]]; then
  : > "_empty.bed"          # create empty output file
  exit 0
fi

# Merge intervals within gap_BP to produce file for chunking
bedtools merge -i combined_counts.bed -d "$GAP_BP" -c 4 -o sum > intervals_with_counts.bed

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

  # merge any abutting intervals
  out="_${pad}${hash}.bed"
  cut -f1-4 "$i" \
  | bedtools merge -i stdin > "$out"
  
  # report size
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$i")
  counts=$(awk '{s+=$4} END{print s+0}' "$i")
  echo "${out}: ${bases} genomic bases, ${counts} summed counts"

  rm $i
done
