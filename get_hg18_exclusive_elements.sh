#!/bin/bash

## hard coded script for qsubbing purposes
cd ~/data/refgenomes/hg18/hg18_refseq/
bedtools complement -i parsed_promoters.merged.bed -g ../hg18.chrom.sizes > non_promoters.bed
bedtools intersect -a parsed_exons.merged.bed -b non_promoters.bed > np_exons.bed

bedtools complement -i parsed_exons.merged.bed -g ../hg18.chrom.sizes > non_exons.bed
bedtools intersect -a parsed_introns.merged.bed -b non_promoters.bed > np_introns.bed
bedtools intersect -a np_introns.bed -b non_exons.bed > npe_introns.bed

bedtools complement -i parsed_introns.merged.bed -g ../hg18.chrom.sizes > non_introns.bed
bedtools intersect -a parsed_repeats.merged.bed -b non_promoters.bed > np_repeats.bed
bedtools intersect -a np_repeats.bed -b non_exons.bed > npe_repeats.bed
bedtools intersect -a npe_repeats.bed -b non_introns.bed > npei_repeats.bed
