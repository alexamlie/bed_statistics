#!/bin/bash

## generate_full_utr_exclusive_elements.sh
## alex amlie-wolf 12-04-2015
## a script to get exclusive elements including 5' and 3' UTRs

## requires bedtools

## currently using the hierarchy:
## 5' UTR exon > 5' UTR intron > 3' UTR exon > 3' UTR intron > promoter > exon > intron > repeat
if [ $# == 3 ]; then
    CHROM_SIZES=$1
    INFOLDER=$2
    OUTFOLDER=$3

    ## assumes that the input folder is organized the same way as the
    ## 12_04_15_protein_coding_file_generation.sh script
    mkdir -p ${OUTFOLDER}/pos_files/final_files/ ${OUTFOLDER}/neg_files/final_files/

    for STRAND in neg pos; do
	echo "Computing exclusive elements on ${STRAND} strand"
	
	cd ${OUTFOLDER}/${STRAND}_files/
	## 5' UTR exon > 5' UTR intron
	## first copy the input file into the final output folder, to be consistent
	cp ${INFOLDER}/parsed_UTR5_exon.${STRAND}.merged.bed final_files/parsed_${STRAND}_5utr_exons.merged.bed
	bedtools complement -i ${INFOLDER}/parsed_UTR5_exon.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_5utr_exons.bed
	bedtools intersect -a ${INFOLDER}/parsed_UTR5_intron.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > final_files/${STRAND}_n5e_5utr_introns.bed

	## 5' UTR intron > 3' UTR exon
	bedtools complement -i ${INFOLDER}/parsed_UTR5_intron.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_5utr_introns.bed
	bedtools intersect -a ${INFOLDER}/parsed_UTR3_exon.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_3utr_exons.bed
	bedtools intersect -a ${STRAND}_n5e_3utr_exons.bed -b ${STRAND}_non_5utr_introns.bed > final_files/${STRAND}_n5e5i_3utr_exons.bed

	## 3' UTR exon > 3' UTR intron
	bedtools complement -i ${INFOLDER}/parsed_UTR3_exon.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_3utr_exons.bed
	bedtools intersect -a ${INFOLDER}/parsed_UTR3_intron.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_3utr_introns.bed
	bedtools intersect -a ${STRAND}_n5e_3utr_introns.bed -b ${STRAND}_non_5utr_introns.bed > ${STRAND}_n5e5i_3utr_introns.bed
	bedtools intersect -a ${STRAND}_n5e5i_3utr_introns.bed -b ${STRAND}_non_3utr_exons.bed > final_files/${STRAND}_n5e5i3e_3utr_introns.bed

	## 3' UTR intron > promoter
	bedtools complement -i ${INFOLDER}/parsed_UTR3_intron.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_3utr_introns.bed
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_promoters.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_promoters.bed
	bedtools intersect -a ${STRAND}_n5e_promoters.bed -b ${STRAND}_non_5utr_introns.bed > ${STRAND}_n5e5i_promoters.bed
	bedtools intersect -a ${STRAND}_n5e5i_promoters.bed -b ${STRAND}_non_3utr_exons.bed > ${STRAND}_n5e5i3e_promoters.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_promoters.bed -b ${STRAND}_non_3utr_introns.bed > final_files/${STRAND}_n5e5i3e3i_promoters.bed

	## promoter > exon
	bedtools complement -i ${INFOLDER}/parsed_mRNA_promoters.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_promoters.bed
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_exon.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_exons.bed
	bedtools intersect -a ${STRAND}_n5e_exons.bed -b ${STRAND}_non_5utr_introns.bed > ${STRAND}_n5e5i_exons.bed
	bedtools intersect -a ${STRAND}_n5e5i_exons.bed -b ${STRAND}_non_3utr_exons.bed > ${STRAND}_n5e5i3e_exons.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_exons.bed -b ${STRAND}_non_3utr_introns.bed > ${STRAND}_n5e5i3e3i_exons.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_exons.bed -b ${STRAND}_non_promoters.bed > final_files/${STRAND}_n5e5i3e3ip_exons.bed

	## exon > intron
	bedtools complement -i ${INFOLDER}/parsed_mRNA_exon.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_exons.bed
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_intron.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_introns.bed
	bedtools intersect -a ${STRAND}_n5e_introns.bed -b ${STRAND}_non_5utr_introns.bed > ${STRAND}_n5e5i_introns.bed
	bedtools intersect -a ${STRAND}_n5e5i_introns.bed -b ${STRAND}_non_3utr_exons.bed > ${STRAND}_n5e5i3e_introns.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_introns.bed -b ${STRAND}_non_3utr_introns.bed > ${STRAND}_n5e5i3e3i_introns.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_introns.bed -b ${STRAND}_non_promoters.bed > ${STRAND}_n5e5i3e3ip_introns.bed
	bedtools intersect -a ${STRAND}_n5e5i3e3ip_introns.bed -b ${STRAND}_non_exons.bed > final_files/${STRAND}_n5e5i3e3ipe_introns.bed

	## intron > repeat
	bedtools complement -i ${INFOLDER}/parsed_mRNA_intron.${STRAND}.merged.bed -g ${CHROM_SIZES} > ${STRAND}_non_introns.bed
	bedtools intersect -a ${INFOLDER}/parsed_repeats.${STRAND}.merged.bed -b ${STRAND}_non_5utr_exons.bed > ${STRAND}_n5e_repeats.bed
	bedtools intersect -a ${STRAND}_n5e_repeats.bed -b ${STRAND}_non_5utr_introns.bed > ${STRAND}_n5e5i_repeats.bed
	bedtools intersect -a ${STRAND}_n5e5i_repeats.bed -b ${STRAND}_non_3utr_exons.bed > ${STRAND}_n5e5i3e_repeats.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_repeats.bed -b ${STRAND}_non_3utr_introns.bed > ${STRAND}_n5e5i3e3i_repeats.bed
	bedtools intersect -a ${STRAND}_n5e5i3e_repeats.bed -b ${STRAND}_non_promoters.bed > ${STRAND}_n5e5i3e3ip_repeats.bed
	bedtools intersect -a ${STRAND}_n5e5i3e3ip_repeats.bed -b ${STRAND}_non_exons.bed > ${STRAND}_n5e5i3e3ipe_repeats.bed
	bedtools intersect -a ${STRAND}_n5e5i3e3ipe_repeats.bed -b ${STRAND}_non_introns.bed > final_files/${STRAND}_n5e5i3e3ipei_repeats.bed

	echo "Sorting output files"

	for F in final_files/*.bed; do
	    echo $F
	    sort -k1,1V -k2,2n $F > $F.sorted
	    mv $F.sorted $F
	done
	
	cd -
	
    done
else
    echo "Usage: $0 <chromosome sizes file> <input folder> <output folder>"
fi
