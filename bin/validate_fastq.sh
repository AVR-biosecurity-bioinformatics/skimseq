#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = library name
# $4 = file1
# $5 = file2

# Returns 0 (true) if gzip is valid AND does NOT have trailing garbage
check_trailing() {
  local out
  if out=$(gzip -t -- "$1" 2>&1); then
    ! grep -qi 'trailing garbage' <<<"$out"
  else
    return 1
  fi
}

okF=false; okR=false
check_trailing "${4}" && okF=true
check_trailing "${5}" && okR=true

# Ensure number of lines are equal
NL_F=$(zcat "${4}" | wc -l | awk '{print $1}' )
NL_R=$(zcat "${5}" | wc -l | awk '{print $1}' )

if $okF && $okR && [[ "$NL_F" -eq "$NL_R" ]]; then
  STATUS=PASS
else
  STATUS=FAIL
fi

# Extract lane, fcid, and platform information from header of first read for read group setup
READ_HEADER=$(zcat ${4} | head -n 1 | sed 's#/1$##' )

# Check if its SRA format data - which doesnt contain FCID and LANE
if [[ $READ_HEADER == @SRR* ]]; then
    # SRA data - Use placeholder FCID and LANE
    FCID=SRA
    LANE=SRA
else
    FCID=$(echo ${READ_HEADER} | cut -d ':' -f 3) #Read flow cell ID
    LANE=$(echo ${READ_HEADER} | cut -d ':' -f 4) #Read lane number 
fi

# TODO: Automatically detect this
# should use "DNBSEQ (MGI/BGI)" for MGI
PLATFORM=ILLUMINA 

# Print information to STDOUT to be split into tuples using nextflow logic
echo -e "${FCID}\t${LANE}\t${PLATFORM}\t${STATUS}"
