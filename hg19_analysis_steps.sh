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







