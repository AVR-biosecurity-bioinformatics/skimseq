#!/bin/bash
set -euo pipefail
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = ref_genome
# $4 = counts_per_chunk
# $5 = split_overweight
# $6 = min_interval_gap
# $7 = include_bed
# $8 = include_zero

# This script uses bed files with or without a counts column to assign bed intervals to seperate chunks for parallel processing
# For Haplotypecaller, input is a single bed with per-base counts in column 4
# For GenotypeGVCFS, input is multiple bed files with no counts
# For Mpileup, input is multiple bed files with per-base counts in column 4

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
  : > "_empty.bed.gz"          # create empty output file
  : > "_empty.bed.gz.tbi"          # create empty output file
  exit 0
fi

btmp="${TMPDIR}/tmp.bed"

SORTED_DIR="${TMPDIR}/sorted"
BEDOPS_INPUT="${TMPDIR}/input.bed5"
MERGED_EXPANDED="${TMPDIR}/merged.expanded.bed"
MERGED_SUMS="${TMPDIR}/merged.sums.bed"

mkdir -p "$SORTED_DIR"
export SORTED_DIR

# Extract included regions, convert BED4 counts to BED5 scores, and sort.
xargs \
    -r \
    -P "${CPUS}" \
    -n 1 \
    bash -c '
        set -euo pipefail

        f="$1"
        base=$(basename "$f" .bed.gz)
        out="${SORTED_DIR}/${base}.sorted.bed5"

        tabix "$f" -R contigs.bed \
        | awk -F "\t" -v OFS="\t" -v file="$f" '"'"'
            NF < 4 {
                printf "ERROR: %s line %d has fewer than four columns: %s\n",
                       file, NR, $0 > "/dev/stderr"
                exit 1
            }

            $1 == "" ||
            $2 !~ /^[0-9]+$/ ||
            $3 !~ /^[0-9]+$/ ||
            $3 <= $2 {
                printf "ERROR: %s line %d has invalid BED coordinates: %s\n",
                       file, NR, $0 > "/dev/stderr"
                exit 1
            }

            $4 !~ /^[-+]?[0-9]*([.][0-9]+)?([eE][-+]?[0-9]+)?$/ {
                printf "ERROR: %s line %d has a non-numeric count: %s\n",
                       file, NR, $0 > "/dev/stderr"
                exit 1
            }

            {
                print $1, $2, $3, ".", $4
            }
        '"'"' \
        | LC_ALL=C sort -t $'"'"'\t'"'"' -k1,1 -k2,2n -k3,3n \
        > "$out"
    ' _ \
    < counts_files.list

# Ensure at least one sorted BED5 was created.
shopt -s nullglob
sorted_beds=("$SORTED_DIR"/*.sorted.bed5)
shopt -u nullglob

if (( ${#sorted_beds[@]} == 0 )); then
    : > "_empty.bed.gz"
    : > "_empty.bed.gz.tbi"
    exit 0
fi

# Combine all sorted BED5 records once.
bedops --everything "${sorted_beds[@]}" \
    > "$BEDOPS_INPUT"

# Extend interval ends by GAP_BP and merge intervals within that distance.
bedops \
    --range "0:${GAP_BP}" \
    --merge \
    "$BEDOPS_INPUT" \
    > "$MERGED_EXPANDED"

# Sum BED5 scores for each merged interval, then remove the temporary
# right-side extension and format integral sums as integers.
bedmap \
    --echo \
    --sum \
    --delim $'\t' \
    "$MERGED_EXPANDED" \
    "$BEDOPS_INPUT" \
| awk -v gap="${GAP_BP}" -v OFS='\t' '
    ($3 - gap > $2) {
        print $1, $2, $3 - gap, int($4)
    }
' > "$MERGED_SUMS"

# Restore the contig order defined by the reference index.
bedtools sort \
    -i "$MERGED_SUMS" \
    -g "${3}.fai" \
    > "$btmp"

# If ALL_BASES=false: keep only intervals that have any counts (doesnt return whole contigs)
# If ALL_BASES=true: map counts back to the full contig intervals (returns whole contigs)
if [[ "$ALL_BASES" == "false" ]]; then
  mv "$btmp" intervals_with_counts.bed
else
  bedtools map -sorted -g ${3}.fai -a contigs.bed -b "$btmp" -c 4 -o sum -null 0 \
  | awk -v OFS='\t' '$4 != 0' \
  > intervals_with_counts.bed
fi

# Exit early if intervals_with_counts.bed is empty
if [[ ! -s intervals_with_counts.bed ]]; then
  : > "_empty.bed.gz"          # create empty output file
  : > "_empty.bed.gz.tbi"      # create empty output file
  exit 0
fi


# Calculate total counts and number of intervals
TOTAL_COUNTS=$(awk '{s+=$4} END{printf "%.0f\n", s+0}' intervals_with_counts.bed)
N_INTERVALS=$(wc -l < intervals_with_counts.bed)

# Decide number of file splits.
# Do not cap by N_INTERVALS when SPLIT_OVERWEIGHT=true, because intervals can be
# split across multiple chunks during assignment.
if [[ "$TOTAL_COUNTS" -le 0 ]]; then
  K=1
else
  K=$(( (TOTAL_COUNTS + TARGET_COUNTS_PER_CHUNK - 1) / TARGET_COUNTS_PER_CHUNK ))
fi

if [[ "$K" -lt 1 ]]; then
  K=1
fi

if [[ "$SPLIT_OVERWEIGHT" != "true" && "$K" -gt "$N_INTERVALS" ]]; then
  K="$N_INTERVALS"
fi

# Assign intervals to chunks.
#
# If SPLIT_OVERWEIGHT=true:
#   intervals are split during assignment whenever they cross the current
#   chunk count boundary. Counts are assigned proportionally/approximately
#   by genomic length, but chunk totals are kept close to the target.
#
# If SPLIT_OVERWEIGHT=false:
#   intervals are kept intact and assigned using cumulative chunk boundaries.
awk -v outdir="$OUTDIR" \
    -v tot="$TOTAL_COUNTS" \
    -v K="$K" \
    -v split_overweight="$SPLIT_OVERWEIGHT" \
    -v OFS="\t" '
function fname(i) {
  return sprintf("%s/chunk_%d.bed", outdir, i)
}

function emit(c, chr, s, e, wt) {
  if (e > s) {
    printf "%s\t%d\t%d\t%.6f\n", chr, s, e, wt >> fname(c)
  }
}

function next_chunk() {
  close(fname(chunk))
  chunk++
  next_cut = chunk * target
  chunk_sum = 0
}

BEGIN {
  if (K < 1) {
    print "ERROR: K must be >= 1" > "/dev/stderr"
    exit 2
  }

  target = tot / K
  chunk = 1
  cumulative = 0
  chunk_sum = 0
  next_cut = target

  # Create/truncate expected output files up front.
  for (i = 1; i <= K; i++) {
    f = fname(i)
    printf "" > f
    close(f)
  }
}

{
  chr = $1
  start = $2 + 0
  end = $3 + 0
  w = $4 + 0

  len = end - start

  if (len <= 0) {
    next
  }

  # Zero-count intervals do not affect balancing.
  if (w <= 0) {
    emit(chunk, chr, start, end, 0)
    next
  }

  # ------------------------------------------------------------------
  # Mode 1: do not split intervals. Assign by cumulative boundaries.
  # ------------------------------------------------------------------
  if (split_overweight != "true") {
    if (chunk < K && chunk_sum > 0) {
      current_dist = cumulative - next_cut
      next_dist = cumulative + w - next_cut

      if ((current_dist < 0 ? -current_dist : current_dist) <= \
          (next_dist < 0 ? -next_dist : next_dist)) {
        next_chunk()
      }
    }

    emit(chunk, chr, start, end, w)
    cumulative += w
    chunk_sum += w
    next
  }

  # ------------------------------------------------------------------
  # Mode 2: split intervals during assignment.
  # ------------------------------------------------------------------
  pos = start
  rem_w = w
  rem_len = len

  while (rem_w > 1e-9 && pos < end) {

    # Last chunk gets all remaining sequence/counts.
    if (chunk >= K) {
      emit(chunk, chr, pos, end, rem_w)
      cumulative += rem_w
      chunk_sum += rem_w
      rem_w = 0
      pos = end
      break
    }

    budget = next_cut - cumulative

    # Floating point guard for exact/near-exact boundaries.
    if (budget <= 1e-9) {
      next_chunk()
      continue
    }

    # Remaining interval fits in the current chunk.
    if (rem_w <= budget + 1e-9) {
      emit(chunk, chr, pos, end, rem_w)
      cumulative += rem_w
      chunk_sum += rem_w
      rem_w = 0
      pos = end
      break
    }

    # Remaining interval crosses the chunk boundary, so split it.
    frac = budget / rem_w
    seg_len = int(rem_len * frac + 0.5)

    # BED intervals must be integer and non-empty.
    if (seg_len < 1) {
      seg_len = 1
    }

    # Avoid creating a zero-length remainder.
    if (seg_len >= rem_len) {
      emit(chunk, chr, pos, end, rem_w)
      cumulative += rem_w
      chunk_sum += rem_w
      rem_w = 0
      pos = end
      break
    }

    split_pos = pos + seg_len

    # Assign exactly the remaining chunk budget to this segment so that
    # chunk totals remain balanced. This assumes approximate uniform count
    # density across the original interval.
    emit(chunk, chr, pos, split_pos, budget)

    pos = split_pos
    rem_len = end - pos
    rem_w -= budget
    cumulative += budget
    chunk_sum += budget

    next_chunk()
  }
}

END {
  for (i = 1; i <= K; i++) {
    close(fname(i))
  }
}
' intervals_with_counts.bed

# Sanity check - make sure number of bases in input matches output
# Input
#awk '{s+=$3-$2} END{print s}' intervals_with_counts.bed
# Chunks
#for f in chunk_*.bed; do
#  awk '{s+=$3-$2} END{print s}' "$f"
#done | awk '{s+=$1} END{print s}'


# Rename each output file
for i in *chunk_*.bed;do
  # Pad output chunk names
  n=$(basename "$i" | sed -E 's/^chunk_([0-9]+)\.bed/\1/')
  pad=$(printf "%05d" "$n")

  # compute hash of this chunk’s contents
  hash=$(md5sum "$i" | awk '{print $1}')

  # Output just column 1:4 as a gzipped bed
  out="${pad}${hash}.bed.gz"
  cut -f1-4 "$i" | bgzip -c > "$out"
  tabix -p bed "$out"

  # report size
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$i")
  counts=$(awk '{s+=$4} END{print s+0}' "$i")
  echo "${out}: ${bases} genomic bases, ${counts} summed counts"

  rm $i
done
