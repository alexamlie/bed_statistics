#!/bin/sh

## bed_shuffle_seq_stats.sh
## alex amlie-wolf
## takes in a bed file, shuffles that file and computes sequence statistics, then outputs
## assumes hg19 for now

GENOME_SIZES=~/data/refgenomes/hg19/hg19.chrom.sizes
GENOME_FA=~/data/refgenomes/hg19/hg19.fa

if [ $# == 2 ];
then
    INFILE=$1
    OUTFILE=$2
    bedtools shuffle -chrom -i $INFILE -g $GENOME_SIZES | bedtools nuc -fi $GENOME_FA -bed stdin > $2
else
    echo "Usage: $0 INPUT_BED OUTPUT_FILE"
fi
    
