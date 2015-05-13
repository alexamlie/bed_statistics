#!/bin/bash

## bed_seq_content.sh
## alex amlie-Wolf
## april/may 2015
## takes a bed file, computes various sequence statistics. basically a wrapper for bedtools nuc
## two arguments: input and output

FA_FILE=/home/alexaml/data/refgenomes/hg19/hg19.fa
CHROMBED=/home/alexaml/data/refgenomes/hg19/hg19_chromosomes.bed

if [ $# == 2 ];
then
    INFILE=$1
    OUTFILE=$2
    bedtools nuc -fi $FA_FILE -bed $INFILE > $OUTFILE
else
    echo "Usage: $0 <input bed file> <output file>"
fi
    
  
