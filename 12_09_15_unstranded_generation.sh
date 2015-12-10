## 12_09_15_unstranded_generation.sh
## alex amlie-wolf
## recording the steps I'm taking (on scisub) to generate bed files for genomic partition analysis

cd ~/data/refgenomes/hg19/unstranded_partitions/input_files/
## make directories to store the bed files, and also the parsed (and sorted) versions
mkdir -p ../bed_files/ ../parsed_files/

for F in hsa19.splitByType.*; do
    echo $F
    TYPE=`echo $F | cut -d. -f3`
    awk -F $'\t' '{ printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, ($4-1), $5, $10, $7 }' $F > ../bed_files/${TYPE}.bed
done

cd ../bed_files/

for F in *; do
    echo $F    
    ~/code/bed_statistics/parse_bedfile_chrs.sh $F ../parsed_files/parsed_${F}
done

cd ../parsed_files/

## now merge these files together
mkdir -p ../merged_files/

## merge all of them
python ~/code/bed_statistics/merge_bed_files.py *.bed
## now move the merged files
mv *merged.bed ../merged_files/

cd ../merged_files/

for F in *.bed; do
    echo $F
    sort -k1,1V -k2,2n $F > $F.sorted
    mv $F.sorted $F
done

## finally, generate the exclusive elements (need bedtools)
module load bedtools2

~/code/bed_statistics/generate_unstranded_elements.sh \
    ~/data/refgenomes/hg19/hg19.chrom.sizes \
    ~/data/refgenomes/hg19/unstranded_partitions/merged_files/ \
    ~/data/refgenomes/hg19/unstranded_partitions/partition/
