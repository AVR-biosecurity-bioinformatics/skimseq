#!/bin/bash
set -euo pipefail

## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = interval_hash
# $4 = .stderr.log files

# Header with the two new columns first
printf "interval_hash\tchrom\tpos\telapsed_min\tvariants\tvariants_per_min\n" > progress_summary.tsv

# Stream: grep only ProgressMeter lines, awk extracts data rows and prefixes sample/hash
grep -F 'ProgressMeter -' "${5}" \
    | awk -v h="${3}" '
    # match data rows like:
    # ... ProgressMeter -   CM028320.1:49552   0.0   24   534.1
    match($0, /ProgressMeter -[[:space:]]*([[:alnum:]_.]+):([0-9]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9.]+)/, m) {
        printf("%s\t%s\t%s\t%s\t%.3f\t%d\t%.1f\n", s, h, m[1], m[2], m[3], m[4], m[5])
    } ' >> progress_summary.tsv
