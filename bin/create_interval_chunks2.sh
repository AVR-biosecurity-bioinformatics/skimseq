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

# Create contig list in reference order with lengths, including just those in the included intervals file
awk 'NR==FNR{c[$1]=1; next} c[$1]' ${7} ${3}.fai \
| cut -f1,2 > contigs.tsv

# Set up memory and cpu control for parralel operations
TOTAL_CPUS="${1}"
TOTAL_MEM_GB="${2}"

JOBS="$TOTAL_CPUS"

# reserve 25% headroom
SORT_MEM_GB=$(( (TOTAL_MEM_GB * 75 / 100) / JOBS ))
(( SORT_MEM_GB < 1 )) && SORT_MEM_GB=1

# Make a list of counts files 
ls -1 *.counts.bed.gz > counts_files.list

# Function for determining how many bed files (samples) overlap each position of the reference genome, and optionally counting bases
process_contig() {
  chr="$1"
  len="$2"

  out="per_contig/${chr}.bed"
  btmp="${TMPDIR}/${chr}.B.tmp"
  all_bed="${TMPDIR}/${chr}.all_intervals.bed"

  genome_line=$(printf "%s\t%s\n" "$chr" "$len")

  # Concatenate all counts beds for that chromosome and merge by gap_bp within file, then within the joint file
  while read -r f; do
    tabix "$f" "$chr" 2>/dev/null \
    | bedtools merge -i - -d "$GAP_BP" -c 4 -o sum
  done < counts_files.list \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -S "$SORT_MEM" -T "$TMPDIR" \
  | bedtools merge -i - -d "$GAP_BP" -c 4 -o sum > "$btmp"

  # If no data, 
  if [ ! -s "$btmp" ]; then
    : > "$out"
    rm -f "$btmp"
    return 0
  fi

  # Create new bedfile which contains all bases with bed records
  if [[ "$ALL_BASES" == "false" ]]; then
    bedtools genomecov -bg -i "$btmp" -g <(printf "%s" "$genome_line") \
      | cut -f1-3 \
      | bedtools merge -i - -d "$GAP_BP"> "$all_bed"
  else 
    # if all bases is set (for short contigs) count across entire contig
    printf "%s\t0\t%s\n" "$chr" "$len" > "$all_bed"
  fi

  # Map column 4 (counts column) back to data
  bedtools map -sorted -a "$all_bed" -b "$btmp" -g <(printf "%s" "$genome_line") -c 4 -o sum -null 0 > "$out"

  # Clean up
  rm -f "$btmp" "$all_bed"
}

# Export variables for GNU parallel
mkdir -p per_contig
export -f process_contig
export GAP_BP
export MODE
export ALL_BASES
export SORT_MEM="${SORT_MEM_GB}G"
export TMPDIR

# Run per-contig process
parallel --colsep '\t' --jobs "${1}" process_contig {1} {2} :::: contigs.tsv

# Merge per-contig counts
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
  ' combined_counts.bed > intervals_split.bed
else
  cat combined_counts.bed > intervals_split.bed
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
