#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = file1
# $4 = file2

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
check_trailing "$3" && okF=true
check_trailing "$4" && okR=true

# Ensure number of lines are equal
NL_F=$(zcat "$3" | wc -l | awk '{print $1}' )
NL_R=$(zcat "$4" | wc -l | awk '{print $1}' )

if $okF && $okR && [[ "$NL_F" -eq "$NL_R" ]]; then
  STATUS=PASS
else
  STATUS=FAIL
fi

# Print ONLY the status to STDOUT (captured by Nextflow as `stdout`)
echo "$STATUS"