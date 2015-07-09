#!/bin/bash

## hard coded script for qsubbing purposes
## unlike get_hg19_exclusive_elements.sh, this one uses different sets of files for each strand
## see hg19_analysis_steps.sh for details
cd ~/data/refgenomes/hg19/hg19_refseq/

## positive strand:
cd pos_strand_files/
bedtools complement -i parsed_pos_promoters.merged.bed -g ../../hg19.chrom.sizes > pos_non_promoters.bed
bedtools intersect -a parsed_pos_exons.merged.bed -b pos_non_promoters.bed > pos_np_exons.bed

bedtools complement -i parsed_pos_exons.merged.bed -g ../../hg19.chrom.sizes > pos_non_exons.bed
bedtools intersect -a parsed_pos_introns.merged.bed -b pos_non_promoters.bed > pos_np_introns.bed
bedtools intersect -a pos_np_introns.bed -b pos_non_exons.bed > pos_npe_introns.bed

bedtools complement -i parsed_pos_introns.merged.bed -g ../../hg19.chrom.sizes > pos_non_introns.bed
bedtools intersect -a parsed_pos_repeats.merged.bed -b pos_non_promoters.bed > pos_np_repeats.bed
bedtools intersect -a pos_np_repeats.bed -b pos_non_exons.bed > pos_npe_repeats.bed
bedtools intersect -a pos_npe_repeats.bed -b pos_non_introns.bed > pos_npei_repeats.bed

## negative strand:
cd ../neg_strand_files/
bedtools complement -i parsed_neg_promoters.merged.bed -g ../../hg19.chrom.sizes > neg_non_promoters.bed
bedtools intersect -a parsed_neg_exons.merged.bed -b neg_non_promoters.bed > neg_np_exons.bed

bedtools complement -i parsed_neg_exons.merged.bed -g ../../hg19.chrom.sizes > neg_non_exons.bed
bedtools intersect -a parsed_neg_introns.merged.bed -b neg_non_promoters.bed > neg_np_introns.bed
bedtools intersect -a neg_np_introns.bed -b neg_non_exons.bed > neg_npe_introns.bed

bedtools complement -i parsed_neg_introns.merged.bed -g ../../hg19.chrom.sizes > neg_non_introns.bed
bedtools intersect -a parsed_neg_repeats.merged.bed -b neg_non_promoters.bed > neg_np_repeats.bed
bedtools intersect -a neg_np_repeats.bed -b neg_non_exons.bed > neg_npe_repeats.bed
bedtools intersect -a neg_npe_repeats.bed -b neg_non_introns.bed > neg_npei_repeats.bed
