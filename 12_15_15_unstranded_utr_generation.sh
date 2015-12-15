## 12_15_15_unstranded_utr_generation.sh
## alex amlie-wolf 12/15/15
## recording the steps i'm taking (on scisub) to generate unstranded bed files for genomic partition
## analysis, including UTR exons and introns

cd ~/data/refgenomes/hg19/unstranded_partitions/
## make a directory to store the UTR annotations
mkdir -p utr_annotations/

## since all the UTR files were generated in 12_09_15_unstranded_generation.sh, just call the
## hierarchy partitioning script with the optional argument to use UTR exons and introns
module load bedtools2

~/code/bed_statistics/generate_unstranded_elements.sh \
    UTR_exons_introns \
    ~/data/refgenomes/hg19/hg19.chrom.sizes \
    ~/data/refgenomes/hg19/unstranded_partitions/merged_files/ \
    ~/data/refgenomes/hg19/unstranded_partitions/utr_annotations/

