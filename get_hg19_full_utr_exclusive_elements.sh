#!/bin/bash

## get_hg19_full_utr_exclusive_elements.sh
## alex amlie-wolf 07-16-2015
## hard coded script for qsubbing purposes
## unlike get_hg19_exclusive_elements.sh, this one uses different sets of files for each strand,
## and unlike get_hg19_strand_exclusive_elements.sh, this one takes 5' and 3' UTRs into account
## unlike get_hg19_utr_exclusive_elements.sh, this one uses both UTR exons and introns
## see hg19_analysis_steps.sh for details

## currently using the hierarchy:
## 5' UTR exon > 5' UTR intron > 3' UTR exon > 3' UTR intron > promoter > exon > intron > repeat
cd ~/data/refgenomes/hg19/hg19_utr_refseq/full_utrs/

CHROM_SIZES=~/data/refgenomes/hg19/hg19.chrom.sizes

## positive strand:
cd pos_files/

## because this generates so many damn files, organize it
mkdir -p final_files

## symlink the 5' UTR exon into the final files folder for clarity
ln -s ~/data/refgenomes/hg19/hg19_utr_refseq/full_utrs/pos_files/parsed_pos_5utr_exons.merged.bed final_files/

## 5' UTR exon > 5' UTR intron
bedtools complement -i parsed_pos_5utr_exons.merged.bed -g ${CHROM_SIZES} > pos_non_5utr_exons.bed
bedtools intersect -a parsed_pos_5utr_introns.merged.bed -b pos_non_5utr_exons.bed > final_files/pos_n5e_5utr_introns.bed

## 5' UTR intron > 3' UTR exon
bedtools complement -i parsed_pos_5utr_introns.merged.bed -g ${CHROM_SIZES} > pos_non_5utr_introns.bed
bedtools intersect -a parsed_pos_3utr_exons.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_3utr_exons.bed
bedtools intersect -a pos_n5e_3utr_exons.bed -b pos_non_5utr_introns.bed > final_files/pos_n5e5i_3utr_exons.bed

## 3' UTR exon > 3' UTR intron
bedtools complement -i parsed_pos_3utr_exons.merged.bed -g ${CHROM_SIZES} > pos_non_3utr_exons.bed
bedtools intersect -a parsed_pos_3utr_introns.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_3utr_introns.bed
bedtools intersect -a pos_n5e_3utr_introns.bed -b pos_non_5utr_introns.bed > pos_n5e5i_3utr_introns.bed
bedtools intersect -a pos_n5e5i_3utr_introns.bed -b pos_non_3utr_exons.bed > final_files/pos_n5e5i3e_3utr_introns.bed

## 3' UTR intron > promoter
bedtools complement -i parsed_pos_3utr_introns.merged.bed -g ${CHROM_SIZES} > pos_non_3utr_introns.bed
bedtools intersect -a parsed_pos_promoters.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_promoters.bed
bedtools intersect -a pos_n5e_promoters.bed -b pos_non_5utr_introns.bed > pos_n5e5i_promoters.bed
bedtools intersect -a pos_n5e5i_promoters.bed -b pos_non_3utr_exons.bed > pos_n5e5i3e_promoters.bed
bedtools intersect -a pos_n5e5i3e_promoters.bed -b pos_non_3utr_introns.bed > final_files/pos_n5e5i3e3i_promoters.bed

## promoter > exon
bedtools complement -i parsed_pos_promoters.merged.bed -g ${CHROM_SIZES} > pos_non_promoters.bed
bedtools intersect -a parsed_pos_exons.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_exons.bed
bedtools intersect -a pos_n5e_exons.bed -b pos_non_5utr_introns.bed > pos_n5e5i_exons.bed
bedtools intersect -a pos_n5e5i_exons.bed -b pos_non_3utr_exons.bed > pos_n5e5i3e_exons.bed
bedtools intersect -a pos_n5e5i3e_exons.bed -b pos_non_3utr_introns.bed > pos_n5e5i3e3i_exons.bed
bedtools intersect -a pos_n5e5i3e_exons.bed -b pos_non_promoters.bed > final_files/pos_n5e5i3e3ip_exons.bed

## exon > intron
bedtools complement -i parsed_pos_exons.merged.bed -g ${CHROM_SIZES} > pos_non_exons.bed
bedtools intersect -a parsed_pos_introns.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_introns.bed
bedtools intersect -a pos_n5e_introns.bed -b pos_non_5utr_introns.bed > pos_n5e5i_introns.bed
bedtools intersect -a pos_n5e5i_introns.bed -b pos_non_3utr_exons.bed > pos_n5e5i3e_introns.bed
bedtools intersect -a pos_n5e5i3e_introns.bed -b pos_non_3utr_introns.bed > pos_n5e5i3e3i_introns.bed
bedtools intersect -a pos_n5e5i3e_introns.bed -b pos_non_promoters.bed > pos_n5e5i3e3ip_introns.bed
bedtools intersect -a pos_n5e5i3e3ip_introns.bed -b pos_non_exons.bed > final_files/pos_n5e5i3e3ipe_introns.bed

## intron > repeat
bedtools complement -i parsed_pos_introns.merged.bed -g ${CHROM_SIZES} > pos_non_introns.bed
bedtools intersect -a parsed_pos_repeats.merged.bed -b pos_non_5utr_exons.bed > pos_n5e_repeats.bed
bedtools intersect -a pos_n5e_repeats.bed -b pos_non_5utr_introns.bed > pos_n5e5i_repeats.bed
bedtools intersect -a pos_n5e5i_repeats.bed -b pos_non_3utr_exons.bed > pos_n5e5i3e_repeats.bed
bedtools intersect -a pos_n5e5i3e_repeats.bed -b pos_non_3utr_introns.bed > pos_n5e5i3e3i_repeats.bed
bedtools intersect -a pos_n5e5i3e_repeats.bed -b pos_non_promoters.bed > pos_n5e5i3e3ip_repeats.bed
bedtools intersect -a pos_n5e5i3e3ip_repeats.bed -b pos_non_exons.bed > pos_n5e5i3e3ipe_repeats.bed
bedtools intersect -a pos_n5e5i3e3ipe_repeats.bed -b pos_non_introns.bed > final_files/pos_n5e5i3e3ipei_repeats.bed

## ----------------------------------------------
## negative strand:
cd ../neg_files/

mkdir -p final_files

## symlink the 5' UTR exon into the final files folder for clarity
ln -s ~/data/refgenomes/hg19/hg19_utr_refseq/full_utrs/neg_files/parsed_neg_5utr_exons.merged.bed final_files/

## 5' UTR exon > 5' UTR intron
bedtools complement -i parsed_neg_5utr_exons.merged.bed -g ${CHROM_SIZES} > neg_non_5utr_exons.bed
bedtools intersect -a parsed_neg_5utr_introns.merged.bed -b neg_non_5utr_exons.bed > final_files/neg_n5e_5utr_introns.bed

## 5' UTR intron > 3' UTR exon
bedtools complement -i parsed_neg_5utr_introns.merged.bed -g ${CHROM_SIZES} > neg_non_5utr_introns.bed
bedtools intersect -a parsed_neg_3utr_exons.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_3utr_exons.bed
bedtools intersect -a neg_n5e_3utr_exons.bed -b neg_non_5utr_introns.bed > final_files/neg_n5e5i_3utr_exons.bed

## 3' UTR exon > 3' UTR intron
bedtools complement -i parsed_neg_3utr_exons.merged.bed -g ${CHROM_SIZES} > neg_non_3utr_exons.bed
bedtools intersect -a parsed_neg_3utr_introns.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_3utr_introns.bed
bedtools intersect -a neg_n5e_3utr_introns.bed -b neg_non_5utr_introns.bed > neg_n5e5i_3utr_introns.bed
bedtools intersect -a neg_n5e5i_3utr_introns.bed -b neg_non_3utr_exons.bed > final_files/neg_n5e5i3e_3utr_introns.bed

## 3' UTR intron > promoter
bedtools complement -i parsed_neg_3utr_introns.merged.bed -g ${CHROM_SIZES} > neg_non_3utr_introns.bed
bedtools intersect -a parsed_neg_promoters.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_promoters.bed
bedtools intersect -a neg_n5e_promoters.bed -b neg_non_5utr_introns.bed > neg_n5e5i_promoters.bed
bedtools intersect -a neg_n5e5i_promoters.bed -b neg_non_3utr_exons.bed > neg_n5e5i3e_promoters.bed
bedtools intersect -a neg_n5e5i3e_promoters.bed -b neg_non_3utr_introns.bed > final_files/neg_n5e5i3e3i_promoters.bed

## promoter > exon
bedtools complement -i parsed_neg_promoters.merged.bed -g ${CHROM_SIZES} > neg_non_promoters.bed
bedtools intersect -a parsed_neg_exons.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_exons.bed
bedtools intersect -a neg_n5e_exons.bed -b neg_non_5utr_introns.bed > neg_n5e5i_exons.bed
bedtools intersect -a neg_n5e5i_exons.bed -b neg_non_3utr_exons.bed > neg_n5e5i3e_exons.bed
bedtools intersect -a neg_n5e5i3e_exons.bed -b neg_non_3utr_introns.bed > neg_n5e5i3e3i_exons.bed
bedtools intersect -a neg_n5e5i3e_exons.bed -b neg_non_promoters.bed > final_files/neg_n5e5i3e3ip_exons.bed

## exon > intron
bedtools complement -i parsed_neg_exons.merged.bed -g ${CHROM_SIZES} > neg_non_exons.bed
bedtools intersect -a parsed_neg_introns.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_introns.bed
bedtools intersect -a neg_n5e_introns.bed -b neg_non_5utr_introns.bed > neg_n5e5i_introns.bed
bedtools intersect -a neg_n5e5i_introns.bed -b neg_non_3utr_exons.bed > neg_n5e5i3e_introns.bed
bedtools intersect -a neg_n5e5i3e_introns.bed -b neg_non_3utr_introns.bed > neg_n5e5i3e3i_introns.bed
bedtools intersect -a neg_n5e5i3e_introns.bed -b neg_non_promoters.bed > neg_n5e5i3e3ip_introns.bed
bedtools intersect -a neg_n5e5i3e3ip_introns.bed -b neg_non_exons.bed > final_files/neg_n5e5i3e3ipe_introns.bed

## intron > repeat
bedtools complement -i parsed_neg_introns.merged.bed -g ${CHROM_SIZES} > neg_non_introns.bed
bedtools intersect -a parsed_neg_repeats.merged.bed -b neg_non_5utr_exons.bed > neg_n5e_repeats.bed
bedtools intersect -a neg_n5e_repeats.bed -b neg_non_5utr_introns.bed > neg_n5e5i_repeats.bed
bedtools intersect -a neg_n5e5i_repeats.bed -b neg_non_3utr_exons.bed > neg_n5e5i3e_repeats.bed
bedtools intersect -a neg_n5e5i3e_repeats.bed -b neg_non_3utr_introns.bed > neg_n5e5i3e3i_repeats.bed
bedtools intersect -a neg_n5e5i3e_repeats.bed -b neg_non_promoters.bed > neg_n5e5i3e3ip_repeats.bed
bedtools intersect -a neg_n5e5i3e3ip_repeats.bed -b neg_non_exons.bed > neg_n5e5i3e3ipe_repeats.bed
bedtools intersect -a neg_n5e5i3e3ipe_repeats.bed -b neg_non_introns.bed > final_files/neg_n5e5i3e3ipei_repeats.bed
