#!/bin/bash

## parse_split_utr_hierarchy.sh
## alex amlie-wolf 07-16-2015
## a script that takes in the entrywise output from entrywise_bed_coverage.py and reports the number
## and proportion of entries in each type of genomic element, assuming that the coverage script
## was run using the split UTR (exons and introns) analysis

## this script assumes a hierarchy:
## 5' UTR exon > 5' UTR intron > 3' UTR exon > 3' UTR intron > promoter > exon > intron > repeat

if [ $# == 1 ]; then
    INFILE=$1
    ## count the exclusive elements in each class
    ## 5' UTR exons:
    NUM_FP_EXON=`tail -n +2 $INFILE | awk '{if ($4>0) print $0}' | wc -l`
    ## 5' UTR introns:
    NUM_FP_INTRON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6>0) print $0}' | wc -l`
    ## 3' UTR exons:
    NUM_TP_EXON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8>0) print $0}' | wc -l`
    ## 3' UTR introns:
    NUM_TP_INTRON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10>0) print $0}' | wc -l`
    ## promoter:
    NUM_PROMOTER=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0 && $12>0) print $0}' | wc -l`
    ## exon:
    NUM_EXON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0 && $12==0 && $14>0) print $0}' | wc -l`
    ## intron:
    NUM_INTRON=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0 && $12==0 && $14==0 && $16>0) print $0}' | wc -l`
    ## repeat:
    NUM_REPEAT=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0 && $12==0 && $14==0 && $16==00 && $18>0) print $0}' | wc -l`
    ## intergenic
    NUM_INTER=`tail -n +2 $INFILE | awk '{if ($4==0 && $6==0 && $8==0 && $10==0 && $12==0 && $14==0 && $16==00 && $18==00) print $0}' | wc -l`
    ## get the total
    TOTAL_NUM=`tail -n +2 $INFILE | wc -l`

    FP_EXON_PROP=`echo ${NUM_FP_EXON}/${TOTAL_NUM} | bc -l`
    FP_INTRON_PROP=`echo ${NUM_FP_INTRON}/${TOTAL_NUM} | bc -l`
    TP_EXON_PROP=`echo ${NUM_TP_EXON}/${TOTAL_NUM} | bc -l`
    TP_INTRON_PROP=`echo ${NUM_TP_INTRON}/${TOTAL_NUM} | bc -l`    
    PROM_PROP=`echo ${NUM_PROMOTER}/${TOTAL_NUM} | bc -l`
    EXON_PROP=`echo ${NUM_EXON}/${TOTAL_NUM} | bc -l`
    INTRON_PROP=`echo ${NUM_INTRON}/${TOTAL_NUM} | bc -l`
    REPEAT_PROP=`echo ${NUM_REPEAT}/${TOTAL_NUM} | bc -l`
    INTER_PROP=`echo ${NUM_INTER}/${TOTAL_NUM} | bc -l`
    TOTAL_PROP=`echo "(${NUM_FP_EXON}+${NUM_FP_INTRON}+${NUM_TP_EXON}+${NUM_TP_INTRON}+${NUM_PROMOTER}+${NUM_EXON}+${NUM_INTRON}+${NUM_REPEAT}+${NUM_INTER})/${TOTAL_NUM}" | bc -l`
    
    printf "Summary of file %s:\n" "$INFILE"
    printf "Type\tNumber\tProportion\n"
    printf "5' UTR exons\t%d\t%.5f\n" "$NUM_FP_EXON" "$FP_EXON_PROP"
    printf "5' UTR introns\t%d\t%.5f\n" "$NUM_FP_INTRON" "$FP_INTRON_PROP"    
    printf "3' UTR exons\t%d\t%.5f\n" "$NUM_TP_EXON" "$TP_EXON_PROP"
    printf "3' UTR introns\t%d\t%.5f\n" "$NUM_TP_INTRON" "$TP_INTRON_PROP"    
    printf "Promoters\t%d\t%.5f\n" "$NUM_PROMOTER" "$PROM_PROP"
    printf "Exons\t%d\t%.5f\n" "$NUM_EXON" "$EXON_PROP"
    printf "Introns\t%d\t%.5f\n" "$NUM_INTRON" "$INTRON_PROP"
    printf "Repeats\t%d\t%.5f\n" "$NUM_REPEAT" "$REPEAT_PROP"
    printf "Intergenic\t%d\t%.5f\n" "$NUM_INTER" "$INTER_PROP"
    printf "Total\t%d\t%.5f\n" "$TOTAL_NUM" "$TOTAL_PROP"
else
    echo "Usage: $0 INPUT_FILE"
fi
