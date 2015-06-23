#!/bin/bash

## bed_to_fasta.sh
## alex amlie-Wolf
## 06-23-2015
## takes a bed file and a genome fasta file and generates an output fasta
## three arguments: input and output

if [ $# == 3 ];
then
    INFILE=$1
    FA_FILE=$2
    OUTFILE=$3
    bedtools getfasta -name -fi $FA_FILE -bed $INFILE -fo $OUTFILE 
else
    echo "Usage: $0 <input bed file> <output file>"
fi
