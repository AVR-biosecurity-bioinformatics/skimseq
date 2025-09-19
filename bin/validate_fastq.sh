#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = file1
# $4 = file2

# Calculate number of reads in forward and reverse fastqs
N_READS_F=$( seqtk size $3 | cut -f1 )
N_READS_R=$( seqtk size $3 | cut -f1 )

if [[ "$N_READS_F" == "$N_READS_R" ]]; then
    # If number of forward and reverse reads match, status = PASS
    STATUS=PASS
else
    # If number of forward and reverse reads dont match, status = FAIL
    STATUS=FAIL
fi

# Output PASS or FAIL for this sample
echo "${2},${STATUS}" > status.csv