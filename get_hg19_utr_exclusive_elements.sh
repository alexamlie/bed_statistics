#!/bin/bash

## get_hg19_utr_exclusive_elements.sh
## alex amlie-wolf 07-15-2015
## hard coded script for qsubbing purposes
## unlike get_hg19_exclusive_elements.sh, this one uses different sets of files for each strand,
## and unlike get_hg19_strand_exclusive_elements.sh, this one takes 5' and 3' UTRs into account
## see hg19_analysis_steps.sh for details

## currently using the hierarchy:
## 5' UTR > promoter > exon > intron > 3' UTR > repeat
cd ~/data/refgenomes/hg19/hg19_utr_refseq/

CHROM_SIZES=~/data/refgenomes/hg19/hg19.chrom.sizes

## positive strand:
cd pos_files/
## 5' UTR > promoter
bedtools complement -i parsed_pos_5utr.merged.bed -g ${CHROM_SIZES} > pos_non_5utr.bed
bedtools intersect -a parsed_pos_promoters.merged.bed -b pos_non_5utr.bed > pos_n5_promoters.bed

## promoter > exon
bedtools complement -i parsed_pos_promoters.merged.bed -g ${CHROM_SIZES} > pos_non_promoters.bed
bedtools intersect -a parsed_pos_exons.merged.bed -b pos_non_5utr.bed > pos_n5_exons.bed
bedtools intersect -a pos_n5_exons.bed -b pos_non_promoters.bed > pos_n5p_exons.bed

## exon > intron
bedtools complement -i parsed_pos_exons.merged.bed -g ${CHROM_SIZES} > pos_non_exons.bed
bedtools intersect -a parsed_pos_introns.merged.bed -b pos_non_5utr.bed > pos_n5_introns.bed
bedtools intersect -a pos_n5_introns.bed -b pos_non_promoters.bed > pos_n5p_introns.bed
bedtools intersect -a pos_n5p_introns.bed -b pos_non_exons.bed > pos_n5pe_introns.bed

## intron > 3' UTR
bedtools complement -i parsed_pos_introns.merged.bed -g ${CHROM_SIZES} > pos_non_introns.bed
bedtools intersect -a parsed_pos_3utr.merged.bed -b pos_non_5utr.bed > pos_n5_3utr.bed
bedtools intersect -a pos_n5_3utr.bed -b pos_non_promoters.bed > pos_n5p_3utr.bed
bedtools intersect -a pos_n5p_3utr.bed -b pos_non_exons.bed > pos_n5pe_3utr.bed
bedtools intersect -a pos_n5pe_3utr.bed -b pos_non_introns.bed > pos_n5pei_3utr.bed

## 3' UTR > repeat
bedtools complement -i parsed_pos_3utr.merged.bed -g ${CHROM_SIZES} > pos_non_3utr.bed
bedtools intersect -a parsed_pos_repeats.merged.bed -b pos_non_5utr.bed > pos_n5_repeats.bed
bedtools intersect -a pos_n5_repeats.bed -b pos_non_promoters.bed > pos_n5p_repeats.bed
bedtools intersect -a pos_n5p_repeats.bed -b pos_non_exons.bed > pos_n5pe_repeats.bed
bedtools intersect -a pos_n5pe_repeats.bed -b pos_non_introns.bed > pos_n5pei_repeats.bed
bedtools intersect -a pos_n5pei_repeats.bed -b pos_non_3utr.bed > pos_n5pei3_repeats.bed

## negative strand:
cd ../neg_files/
## 5' UTR > promoter
bedtools complement -i parsed_neg_5utr.merged.bed -g ${CHROM_SIZES} > neg_non_5utr.bed
bedtools intersect -a parsed_neg_promoters.merged.bed -b neg_non_5utr.bed > neg_n5_promoters.bed

## promoter > exon
bedtools complement -i parsed_neg_promoters.merged.bed -g ${CHROM_SIZES} > neg_non_promoters.bed
bedtools intersect -a parsed_neg_exons.merged.bed -b neg_non_5utr.bed > neg_n5_exons.bed
bedtools intersect -a neg_n5_exons.bed -b neg_non_promoters.bed > neg_n5p_exons.bed

## exon > intron
bedtools complement -i parsed_neg_exons.merged.bed -g ${CHROM_SIZES} > neg_non_exons.bed
bedtools intersect -a parsed_neg_introns.merged.bed -b neg_non_5utr.bed > neg_n5_introns.bed
bedtools intersect -a neg_n5_introns.bed -b neg_non_promoters.bed > neg_n5p_introns.bed
bedtools intersect -a neg_n5p_introns.bed -b neg_non_exons.bed > neg_n5pe_introns.bed

## intron > 3' UTR
bedtools complement -i parsed_neg_introns.merged.bed -g ${CHROM_SIZES} > neg_non_introns.bed
bedtools intersect -a parsed_neg_3utr.merged.bed -b neg_non_5utr.bed > neg_n5_3utr.bed
bedtools intersect -a neg_n5_3utr.bed -b neg_non_promoters.bed > neg_n5p_3utr.bed
bedtools intersect -a neg_n5p_3utr.bed -b neg_non_exons.bed > neg_n5pe_3utr.bed
bedtools intersect -a neg_n5pe_3utr.bed -b neg_non_introns.bed > neg_n5pei_3utr.bed

## 3' UTR > repeat
bedtools complement -i parsed_neg_3utr.merged.bed -g ${CHROM_SIZES} > neg_non_3utr.bed
bedtools intersect -a parsed_neg_repeats.merged.bed -b neg_non_5utr.bed > neg_n5_repeats.bed
bedtools intersect -a neg_n5_repeats.bed -b neg_non_promoters.bed > neg_n5p_repeats.bed
bedtools intersect -a neg_n5p_repeats.bed -b neg_non_exons.bed > neg_n5pe_repeats.bed
bedtools intersect -a neg_n5pe_repeats.bed -b neg_non_introns.bed > neg_n5pei_repeats.bed
bedtools intersect -a neg_n5pei_repeats.bed -b neg_non_3utr.bed > neg_n5pei3_repeats.bed
