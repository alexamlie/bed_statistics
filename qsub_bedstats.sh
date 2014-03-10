#!/bin/bash

## qsub_bedstats.sh
## Alex Amlie-Wolf
## A simple wrapper script to run my bed statistics 
## for ease, takes in the files you want to use as arguments

if [ $# -lt 3 ]
then
    echo "USAGE: $0 reffile bedfile output"
    exit 1
fi

export CODE=/home/alexaml/code/bed_statistics/

python $CODE/bed_statistics.py $1 $2 $3
