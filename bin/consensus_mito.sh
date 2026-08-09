#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = cram
# $4 = ref_genome fasta
# $5 = mito_genome fasta
# $6 = mito_contig bed
# $7 = numt bed
# $8 = min_vaf
# $9 = min_depth

# Merge real mito and numt beds
cat "${6}" "${7}" \
  | awk 'BEGIN{OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > mt_numt.bed

# Define variables
REF=${4}
MT=${5}
MINVAF=${8}
MINDEPTH=${9}

# Extract all reads aligning mitochondrial contig and numt regions from merged bam and convert to FASTQ 
samtools view -b -T "$REF" --regions-file mt_numt.bed ${3} \
| samtools collate -Ou - \
| samtools fastq -o mito.fq 
    
# Fuction for counting alleles using samtools consensus pileup
# NOTE: this has been tested robustly
count_mpileup_alleles() {
  # mpileup columns: rname pos ref depth bases quals
  # -aa: emit all positions (including depth 0)
  # -Q/-q: base and mapQ filters; set as you prefer
  # -x ignore overlap removal
  # -A count orphans
  samtools mpileup -aa -A -q 0 -Q 0 -B --max-depth 10000 -f "$MT" - \
  | awk -v OFS="\t" -v MINVAF="$MINVAF" -v MINDEPTH="$MINDEPTH" '
    {
      rname = $1
      pos   = $2
      ref   = toupper($3)
      bases = $5

      # mpileup uses: . , for ref matches; ^. for read start; $ for read end
      # and +nSEQ / -nSEQ indel annotations inside the bases string.
      # We must strip control characters & indel payload before counting.

      # If there are no reads, bases may be "*" or empty depending on flags;
      # handle robustly after counting.

      # Remove start-of-read markers and the following MQ char
      gsub(/\^./, "", bases)
      # Remove end-of-read markers
      gsub(/\$/, "", bases)

      # Convert ref matches (.,) into the reference base so they count as A/C/G/T
      if (ref ~ /^[ACGT]$/) {
        gsub(/\./, ref, bases)
        gsub(/,/, ref, bases)
      } else {
        # If ref is N/other, convert .,/ to N so they don’t inflate A/C/G/T
        gsub(/\./, "N", bases)
        gsub(/,/, "N", bases)
      }

      bases = toupper(bases)

      # Optional: count insertion events before stripping them
      # (each +<len><seq> counts as one event)
      ins_events = 0
      tmp2 = bases
      while (match(tmp2, /\+[0-9]+/)) {
        ins_events++
        n = substr(tmp2, RSTART+1, RLENGTH-1) + 0
        tmp2 = substr(tmp2, 1, RSTART-1) substr(tmp2, RSTART+RLENGTH+n)
      }

      # Strip all indel annotations (+nSEQ and -nSEQ) from bases
      # This leaves base calls at this position plus deletion placeholders (*)
      while (match(bases, /[+-][0-9]+/)) {
        n = substr(bases, RSTART+1, RLENGTH-1) + 0
        bases = substr(bases, 1, RSTART-1) substr(bases, RSTART+RLENGTH+n)
      }

      # Now count A/C/G/T and deletion placeholders (*)
      tmp=bases; A = gsub(/A/, "", tmp)
      tmp=bases; C = gsub(/C/, "", tmp)
      tmp=bases; G = gsub(/G/, "", tmp)
      tmp=bases; T = gsub(/T/, "", tmp)
      tmp=bases; del = gsub(/\*/, "", tmp)

      depth = A + C + G + T + del
      nth = 0

      # Handle positions with no usable observations
      if (depth == 0) {
        consensus = "N"
        vaf = 0
        print rname, pos, nth, consensus, depth, A, C, G, T, del #, ins_events
        next
      }

      # Determine major allele among A/C/G/T (ignore del for consensus call)
      consensus="A"; max=A; tie=0
      if (C > max) {consensus="C"; max=C; tie=0} else if (C==max && max>0) tie=1
      if (G > max) {consensus="G"; max=G; tie=0} else if (G==max && max>0) tie=1
      if (T > max) {consensus="T"; max=T; tie=0} else if (T==max && max>0) tie=1

      vaf = max / depth

      # Depth/VAF/tie filtering
      if (depth < MINDEPTH || vaf < MINVAF || tie==1) {
        consensus = "N"
      }

      print rname, pos, nth, consensus, depth, A, C, G, T, del #, ins_events
    }'
}

# Align to original reference and count alleles
bwa-mem2 mem -t "$1" -p -a "$MT" mito.fq \
  | samtools sort -@ "$1" -O BAM -o - - \
  | count_mpileup_alleles \
  > allele_counts.txt

# Create shifted mitochondrial reference
MT_SHIFTED=shifted.fa
SHIFT=5000
MT_LEN=$(awk 'NR==1{print $2}' "${MT}.fai")
MT_CONTIG=$(awk 'NR==1{print $1}' "${MT}.fai")

seqkit restart -i $SHIFT ${MT} > $MT_SHIFTED
samtools faidx $MT_SHIFTED
bwa-mem2 index $MT_SHIFTED

# Shifted
bwa-mem2 mem -t "$1" -p -a "$MT_SHIFTED" mito.fq \
  | samtools sort -@ "$1" -O BAM -o - - \
  | count_mpileup_alleles \
  > allele_counts_shifted.txt

# unshift the allele counts
awk -v OFS="\t" -v SHIFT="$SHIFT" -v L="$MT_LEN" -v CONT="$MT_CONTIG" '
  function mod(a,m) { return ((a % m) + m) % m }  # negative-safe modulo
  {
    pos = $2 + 0

    # SeqKit restart semantics:
    #   SHIFT>0  -> start at base SHIFT
    #   SHIFT<0  -> start at base (L+SHIFT+1)
    se = SHIFT
    if (se < 0) se = L + se + 1

    # Map shifted pos -> original pos (1-based)
    pos2 = mod((pos + se - 2), L) + 1

    $1 = CONT
    $2 = pos2
    print
  }' allele_counts_shifted.txt \
| sort -k2,2n -k3,3n \
> allele_counts_unshifted.txt

# Combine original and unshifted allele counts files, prefering the one with the higher counts
awk -v OFS="\t" '
FNR==NR{
  key=$1 FS $2 FS $3
  best_d[key]=$5+0
  best_line[key]=$0
  next
}
{
  key=$1 FS $2 FS $3
  d=$5+0
  if (!(key in best_d) || d > best_d[key]) {
    best_d[key]=d
    best_line[key]=$0
  }
}
END{
  for (k in best_line) print best_line[k]
}
' allele_counts.txt allele_counts_unshifted.txt \
| sort -k2,2n -k3,3n \
> ${3}.allele_counts.txt

# Extract the consensus column and write to fasta
cut -f4 ${3}.allele_counts.txt \
  | tr -d '\n' \
  | sed -e "1i>${2}" -e 's/.\{60\}/&\n/g' \
  > ${3}.mito.fa

# clean up
rm mito.fq 