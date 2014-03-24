#!/bin/bash
#PBS -o qsub_output.out
#PBS -e qsub_error.err

## qsub_bedstats.sh
## Alex Amlie-Wolf
## A simple wrapper script to run my bed statistics 
## for ease, takes in the files you want to use as arguments


if [ $# -lt 4 ]
then
    echo "USAGE: $0 reffile bedfile output log"
    exit 1
fi

export CODE=/home/alexaml/code/bed_statistics/

# move working directory to desired place
cd `dirname $3`

# use python 2.7
source ~/virtualenvs/py2.7/bin/activate
python $CODE/bed_statistics.py $1 $2 $3 $4
deactivate

# move back (just in case?)
cd -
