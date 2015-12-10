#!/bin/sh

## perform_bed_partition_analysis.sh
## alex amlie-wolf 12-09-2015
## a script that takes in a bed file and performs genomic partition analysis
## split if it's stranded or not

if [ $# == 3 ]; then
    INBED=$1
    OUTDIR=$2
    PARTITION_DIR=$3
    
    mkdir -p ${OUTDIR}

    NCOL=`awk '{print NF}' $INBED | sort -u`
    BED_NAME=`basename ${INBED%.bed}`

    ## sort the bed file:
    sort -k1,1V -k2,2n ${INBED} > ${OUTDIR}/${BED_NAME}_sorted.bed
    
    if [ ${NCOL} -ge 6 ]; then
	awk '{if ($6=="+") print $0}' ${OUTDIR}/${BED_NAME}_sorted.bed | sort -k1,1V -k2,2n \
	    > ${OUTDIR}/pos_${BED_NAME}.bed
	awk '{if ($6=="-") print $0}' ${OUTDIR}/${BED_NAME}_sorted.bed | sort -k1,1V -k2,2n \
	    > ${OUTDIR}/neg_${BED_NAME}.bed
	
	## hard coded code dir
	python ~/code/bed_statistics/entrywise_bed_coverage.py \
	    --full_utrs \
	    ${PARTITION_DIR}/pos_files/final_files/parsed_pos_5utr_exons.merged.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e_5utr_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i_3utr_exons.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e_3utr_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3i_promoters.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ip_exons.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ipe_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ipei_repeats.bed \
	    ${OUTDIR}/pos_${BED_NAME}.bed \
	    ${OUTDIR}/pos_${BED_NAME}_entry.txt ${OUTDIR}/pos_${BED_NAME}_summary.txt
	
	python ~/code/bed_statistics/entrywise_bed_coverage.py \
	    --full_utrs \
	    ${PARTITION_DIR}/neg_files/final_files/parsed_neg_5utr_exons.merged.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e_5utr_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i_3utr_exons.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e_3utr_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3i_promoters.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ip_exons.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ipe_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ipei_repeats.bed \
	    ${OUTDIR}/neg_${BED_NAME}.bed \
	    ${OUTDIR}/neg_${BED_NAME}_entry.txt ${OUTDIR}/neg_${BED_NAME}_summary.txt

	## write the summary to file
	~/code/bed_statistics/parse_split_utr_hierarchy.sh ${OUTDIR}/pos_${BED_NAME}_entry.txt \
	    > ${OUTDIR}/pos_${BED_NAME}_parsed_hierarchy.txt

	~/code/bed_statistics/parse_split_utr_hierarchy.sh ${OUTDIR}/neg_${BED_NAME}_entry.txt \
	    > ${OUTDIR}/neg_${BED_NAME}_parsed_hierarchy.txt
    ## if we aren't stranded
    else
	echo "Performing non-stranded analysis"
	## hard coded code dir
	python ~/code/bed_statistics/entrywise_bed_coverage.py \
	    --full_utrs \
	    ${PARTITION_DIR}/pos_files/final_files/parsed_pos_5utr_exons.merged.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e_5utr_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i_3utr_exons.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e_3utr_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3i_promoters.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ip_exons.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ipe_introns.bed \
	    ${PARTITION_DIR}/pos_files/final_files/pos_n5e5i3e3ipei_repeats.bed \
	    ${OUTDIR}/${BED_NAME}_sorted.bed \
	    ${OUTDIR}/pos_${BED_NAME}_entry.txt ${OUTDIR}/pos_${BED_NAME}_summary.txt
	
	python ~/code/bed_statistics/entrywise_bed_coverage.py \
	    --full_utrs \
	    ${PARTITION_DIR}/neg_files/final_files/parsed_neg_5utr_exons.merged.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e_5utr_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i_3utr_exons.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e_3utr_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3i_promoters.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ip_exons.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ipe_introns.bed \
	    ${PARTITION_DIR}/neg_files/final_files/neg_n5e5i3e3ipei_repeats.bed \
	    ${OUTDIR}/${BED_NAME}_sorted.bed \
	    ${OUTDIR}/neg_${BED_NAME}_entry.txt ${OUTDIR}/neg_${BED_NAME}_summary.txt

	## write the summary to file
	~/code/bed_statistics/parse_split_utr_hierarchy.sh ${OUTDIR}/pos_${BED_NAME}_entry.txt \
	    > ${OUTDIR}/pos_${BED_NAME}_parsed_hierarchy.txt

	~/code/bed_statistics/parse_split_utr_hierarchy.sh ${OUTDIR}/neg_${BED_NAME}_entry.txt \
	    > ${OUTDIR}/neg_${BED_NAME}_parsed_hierarchy.txt
    fi
else
    echo "Usage: $0 <input bed file> <output directory>"
fi
