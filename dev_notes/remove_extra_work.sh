#! /bin/bash

## Modified from https://gist.github.com/photocyte/495848faaba3319c962a575593eaeb55

## Has to be run in directory the nextflow pipeline was run from

## Find work directories essential to the last pipeline run, as absolute paths
nextflow log last > /tmp/preserve_dirs.txt

# Get work directory path
workdir=$(cat /tmp/preserve_dirs.txt | head -1 | sed -E 's|^(.*/work)/.*|\1|')

## Find all work directories, as absolute paths
find "$(readlink -f $workdir)" -maxdepth 2 -type d -path '**/work/*/*' > /tmp/all_dirs.txt

## Concatenate, sort, and count, filtering to those that show up only once (i.e., just once from all_dirs.txt)
grep -Fvxf /tmp/preserve_dirs.txt /tmp/all_dirs.txt > /tmp/to_delete_dirs.txt

## Delete the extraneous work directories
cat /tmp/to_delete_dirs.txt | xargs -r -P 4 -n 1 rm -rf

## Clean up
rm -f /tmp/preserve_dirs.txt /tmp/all_dirs.txt /tmp/to_delete_dirs.txt