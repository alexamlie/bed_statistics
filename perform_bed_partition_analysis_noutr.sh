#!/bin/sh

## perform_bed_partition_analysis_noutr.sh
## alex amlie-wolf 12-09-2015
## a script that takes in a bed file and performs genomic partition analysis
## uses unstranded loci with no UTRs

if [ $# == 3 ]; then
    INBED=$1
    OUTDIR=$2
    PARTITION_DIR=$3
    
    mkdir -p ${OUTDIR}

    BED_NAME=`basename ${INBED%.bed}`

    ## sort the bed file:
    sort -k1,1V -k2,2n ${INBED} > ${OUTDIR}/${BED_NAME}_sorted.bed
    
    ## hard coded code dir
    python ~/code/bed_statistics/entrywise_bed_coverage.py \
    	${PARTITION_DIR}/final_files/parsed_mRNA_promoters.merged.bed \
    	${PARTITION_DIR}/final_files/np_exons.bed \
    	${PARTITION_DIR}/final_files/npe_introns.bed \
    	${PARTITION_DIR}/final_files/npei_repeats.bed \
    	${OUTDIR}/${BED_NAME}_sorted.bed \
    	${OUTDIR}/${BED_NAME}_entry.txt ${OUTDIR}/${BED_NAME}_summary.txt
    
    ## write the summary to file
    ~/code/bed_statistics/parse_entrywise_exclusive_classes.sh ${OUTDIR}/${BED_NAME}_entry.txt \
	> ${OUTDIR}/${BED_NAME}_parsed_hierarchy.txt    
else
    echo "Usage: $0 <input bed file> <output directory> <partition directory>"
fi
