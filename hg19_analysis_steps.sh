## hg19_analysis_steps.sh
## the steps I ran to get promoter, exon, intron, and intergenic files for hg19
## the first step is downloading the 'raw' tables of each element from the UCSC table browser
## the promoters are from refseq genes, 1000bp up, and the exons and introns are refseq
## the repeats are from repeat masker
cd ~/code/bed_statistics/
mkdir -p logs/

## first, parse the input bed files so that they have only chr1-22, X, and Y
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_promoters ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_refseq_1000bpup_promoters.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_promoters.bed
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_exons ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_refseq_exons.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_exons.bed
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_introns ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_refseq_introns.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_introns.bed
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_intergenic ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_repeats.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_repeats.bed

## now merge overlapping entries
qsub -wd ~/code/bed_statistics/logs/ -N merge_hg19_reference ~/code/qsub_wrappers/python_wrapper.sh ~/code/bed_statistics/merge_reference_beds.py ~/data/refgenomes/hg19/hg19_refseq/parsed_promoters.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_exons.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_introns.bed ~/data/refgenomes/hg19/hg19_refseq/parsed_repeats.bed

## now get the exclusive element loci
qsub -wd ~/code/bed_statistics/logs/ -N get_hg19_exclusive_elements ~/code/bed_statistics/get_hg19_exclusive_elements.sh

## -------------
## SPLIT BY STRAND ANALYSIS:
## -------------
cd ~/data/refgenomes/hg19/hg19_refseq/
mkdir -p pos_strand_files/ neg_strand_files/

## take the parsed input files, split into positive and negative
awk '{if ($6=="+") print $0}' parsed_promoters.bed > pos_strand_files/parsed_pos_promoters.bed
awk '{if ($6=="+") print $0}' parsed_exons.bed > pos_strand_files/parsed_pos_exons.bed
awk '{if ($6=="+") print $0}' parsed_introns.bed > pos_strand_files/parsed_pos_introns.bed
awk '{if ($6=="+") print $0}' parsed_repeats.bed > pos_strand_files/parsed_pos_repeats.bed

awk '{if ($6=="-") print $0}' parsed_promoters.bed > neg_strand_files/parsed_neg_promoters.bed
awk '{if ($6=="-") print $0}' parsed_exons.bed > neg_strand_files/parsed_neg_exons.bed
awk '{if ($6=="-") print $0}' parsed_introns.bed > neg_strand_files/parsed_neg_introns.bed
awk '{if ($6=="-") print $0}' parsed_repeats.bed > neg_strand_files/parsed_neg_repeats.bed

## now merge these
qsub -wd ~/code/bed_statistics/logs/ -N merge_hg19_pos_reference ~/code/qsub_wrappers/python_wrapper.sh ~/code/bed_statistics/merge_reference_beds.py ~/data/refgenomes/hg19/hg19_refseq/pos_strand_files/parsed_pos_promoters.bed ~/data/refgenomes/hg19/hg19_refseq/pos_strand_files/parsed_pos_exons.bed ~/data/refgenomes/hg19/hg19_refseq/pos_strand_files/parsed_pos_introns.bed ~/data/refgenomes/hg19/hg19_refseq/pos_strand_files/parsed_pos_repeats.bed

qsub -wd ~/code/bed_statistics/logs/ -N merge_hg19_neg_reference ~/code/qsub_wrappers/python_wrapper.sh ~/code/bed_statistics/merge_reference_beds.py ~/data/refgenomes/hg19/hg19_refseq/neg_strand_files/parsed_neg_promoters.bed ~/data/refgenomes/hg19/hg19_refseq/neg_strand_files/parsed_neg_exons.bed ~/data/refgenomes/hg19/hg19_refseq/neg_strand_files/parsed_neg_introns.bed ~/data/refgenomes/hg19/hg19_refseq/neg_strand_files/parsed_neg_repeats.bed

## finally, generate the exclusive elements for each strand
qsub -wd ~/code/bed_statistics/logs/ -N get_hg19_stranded_exclusive_elements ~/code/bed_statistics/get_hg19_strand_exclusive_elements.sh

## ----------------------------------------------
## UTR analysis (07-15-2015):
## ----------------------------------------------
## we want to add UTR annotations to the partitioning analysis, so here are the steps I used to
## generate the files for this analysis:
cd ~/data/refgenomes/hg19/
mkdir hg19_utr_refseq ## put it in a separate folder
cd hg19_utr_refseq
mkdir -p pos_files/ neg_files/

## link the parsed and merged other files here (strand specific)
for f in ~/data/refgenomes/hg19/hg19_refseq/pos_strand_files/parsed_*.merged.bed; do
    ln -s $f pos_files/
done

for f in ~/data/refgenomes/hg19/hg19_refseq/neg_strand_files/parsed_*.merged.bed; do
    ln -s $f neg_files/
done

## now parse the UTR files:
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_5utr ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_refseq_5utr.bed ~/data/refgenomes/hg19/hg19_utr_refseq/parsed_5utr.bed
qsub -wd ~/code/bed_statistics/logs/ -N hg19_parse_3utr ~/code/bed_statistics/parse_bedfile_chrs.sh ~/data/refgenomes/hg19/hg19_refseq_3utr.bed ~/data/refgenomes/hg19/hg19_utr_refseq/parsed_3utr.bed

## split these by strand, then merge:
awk '{if ($6=="+") print $0}' parsed_5utr.bed > pos_files/parsed_pos_5utr.bed
awk '{if ($6=="+") print $0}' parsed_3utr.bed > pos_files/parsed_pos_3utr.bed
awk '{if ($6=="-") print $0}' parsed_5utr.bed > neg_files/parsed_neg_5utr.bed
awk '{if ($6=="-") print $0}' parsed_5utr.bed > neg_files/parsed_neg_3utr.bed

## merge overlapping entries here
## i have created a new merging script that just takes an arbitrary number of bed files:
qsub -wd ~/code/bed_statistics/logs/ -N merge_hg19_utrs ~/code/qsub_wrappers/python_wrapper.sh ~/code/bed_statistics/merge_bed_files.py ~/data/refgenomes/hg19/hg19_utr_refseq/pos_files/parsed_pos_5utr.bed ~/data/refgenomes/hg19/hg19_utr_refseq/pos_files/parsed_pos_3utr.bed ~/data/refgenomes/hg19/hg19_utr_refseq/neg_files/parsed_neg_5utr.bed ~/data/refgenomes/hg19/hg19_utr_refseq/neg_files/parsed_neg_3utr.bed 

## finally, generate the exclusive elements for each strand (with UTRs)
qsub -wd ~/code/bed_statistics/logs/ -N get_hg19_utr_exclusive_elements ~/code/bed_statistics/get_hg19_utr_exclusive_elements.sh
