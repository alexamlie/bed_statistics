#!/bin/bash

## parse_entrywise_class_hierarchy.sh
## alex amlie-wolf 06-16-15
## a script that takes in the entrywise output from entrywise_bed_coverage.py and reports the number
## and proportion of entries in each type of genomic element. unlike the other parsing scripts,
## this script assumes a hierarchy: promoters > exons > introns > repeat > intergenic, so if there
## is any promoter overlap, the entry is a promoter regardless of exonic and intronic overlap, etc.

if [ $# == 1 ]; then
    INFILE=$1
    ## count the exclusive elements in each class
    ## promoter:
    NUM_PROMOTER=`tail -n +2 $INFILE | awk '{if ($4>0) print $0}' | wc -l`
    ## exon:
    NUM_EXON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6>0) print $0}' | wc -l`
    ## intron:
    NUM_INTRON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8>0) print $0}' | wc -l`
    ## repeat:
    NUM_REPEAT=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10>0) print $0}' | wc -l`
    ## intergenic
    NUM_INTER=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0) print $0}' | wc -l`
    ## get the total
    TOTAL_NUM=`tail -n +2 $INFILE | wc -l`

    PROM_PROP=`echo ${NUM_PROMOTER}/${TOTAL_NUM} | bc -l`
    EXON_PROP=`echo ${NUM_EXON}/${TOTAL_NUM} | bc -l`
    INTRON_PROP=`echo ${NUM_INTRON}/${TOTAL_NUM} | bc -l`
    REPEAT_PROP=`echo ${NUM_REPEAT}/${TOTAL_NUM} | bc -l`
    INTER_PROP=`echo ${NUM_INTER}/${TOTAL_NUM} | bc -l`
    
    printf "Summary of file %s:\n" "$INFILE"
    printf "Type\tNumber\tProportion\n"
    printf "Promoters\t%d\t%.5f\n" "$NUM_PROMOTER" "$PROM_PROP"
    printf "Exons\t%d\t%.5f\n" "$NUM_EXON" "$EXON_PROP"
    printf "Introns\t%d\t%.5f\n" "$NUM_INTRON" "$INTRON_PROP"
    printf "Repeats\t%d\t%.5f\n" "$NUM_REPEAT" "$REPEAT_PROP"
    printf "Intergenic\t%d\t%.5f\n" "$NUM_INTER" "$INTER_PROP"
    printf "Total\t%d\t%.5f\n" "$TOTAL_NUM" "1.0000"
else
    echo "Usage: $0 INPUT_FILE"
fi
