#! /bin/bash
set -euo pipefail

## Modified from https://gist.github.com/photocyte/495848faaba3319c962a575593eaeb55

## Has to be run in directory the nextflow pipeline was run from
tmpdir=$(mktemp -d)
preserve_dirs="${tmpdir}/preserve_dirs.txt"
all_dirs="${tmpdir}/all_dirs.txt"
to_delete_dirs="${tmpdir}/to_delete_dirs.txt"

## Find work directories required by the most recent run
nextflow log last > "$preserve_dirs"

[[ -s "$preserve_dirs" ]] || {
    echo "ERROR: nextflow log last returned no task directories"
    exit 1
}

## Determine work directory
workdir=$(head -n 1 "$preserve_dirs" | sed -E 's|^(.*/work)/.*|\1|')

[[ -n "$workdir" ]] || {
    echo "ERROR: Could not determine work directory"
    exit 1
}

[[ -d "$workdir" ]] || {
    echo "ERROR: Work directory does not exist: $workdir"
    exit 1
}

workdir=$(readlink -f "$workdir")


## Find all Nextflow task directories
find "$workdir" -mindepth 2 -maxdepth 2 -type d > "$all_dirs"

## Determine which directories can be deleted
grep -Fvxf "$preserve_dirs" "$all_dirs" > "$to_delete_dirs" || true

echo "Preserving $(wc -l < "$preserve_dirs") work directories"
echo "Deleting  $(wc -l < "$to_delete_dirs") work directories"

## Delete the extraneous work directories
grep -F "${workdir}/" ${to_delete_dirs} \
    | xargs -r -P 8 -n 1 rm -rf

## Clean up
rm -f ${all_dirs} ${to_delete_dirs} ${preserve_dirs}