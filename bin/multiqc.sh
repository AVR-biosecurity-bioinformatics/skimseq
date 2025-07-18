#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = config file

multiqc \
    --force \
    ${3} \
    .