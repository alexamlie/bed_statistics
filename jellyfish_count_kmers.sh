#!/bin/bash

## jellyfish_count_kmers.sh
## alex amlie-wolf 06-23-15
## uses the jellyfish program to count kmers
## takes in the input fasta file, the desired value of k, and the output file

TMPDIR=~/data/tmp/

if [ $# == 3 ]; then
    INFILE=$1
    KVAL=$2
    OUTFILE=$3

    ## write to a temporary file and then extract
    jellyfish count -m $KVAL -s 100M -t 1 -C -o $TMPDIR/${$}.tmp $INFILE
    ## get the counts
    jellyfish dump -c -t -o $OUTFILE $TMPDIR/${$}.tmp
    jellyfish info $TMPDIR/${$}.tmp > $OUTFILE.info
    rm $TMPDIR/${$}.tmp
else
    echo "Usage: $0 INFILE KVALUE OUTFILE"
fi
