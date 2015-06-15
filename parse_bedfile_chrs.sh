#!/usr/bin/bash

## parse_bedfile_chrs.sh
## alex amlie-wolf 06-11-2015
## takes in bed files, parses them to only have 'normal' chromosomes (1-22, X, and Y) and sorts them
## by chromosome, start position, and strand
## first argument is bed file, second is output

if [ $# == 2 ]; then
    grep "^chr[0-9XY]*[[:space:]]" $1 | awk 'BEGIN{OFS="\t"} {if($3 >= $2) print $0}' | sort -u -t$'\t' -T ~/data/tmp/ -k1,1 -k2,2n -k6,6 > $2
else
    echo "Usage: $0 INPUT_BED OUTPUT_FILE"
fi

