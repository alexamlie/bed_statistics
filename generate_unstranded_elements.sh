if [ $# == 3 ]; then
    CHROM_SIZES=$1
    INFOLDER=$2
    OUTFOLDER=$3

    ## assumes that the input folder is organized the same way as the
    ## 12_04_15_protein_coding_file_generation.sh script
    mkdir -p ${OUTFOLDER}/final_files/

    cd ${OUTFOLDER}

    ## generate all the complements:
    ## promoters > exons > introns > repeats 
    bedtools complement -i ${INFOLDER}/parsed_mRNA_promoters.merged.bed -g ${CHROM_SIZES} > non_promoters.bed
    bedtools complement -i ${INFOLDER}/parsed_mRNA_exon.merged.bed -g ${CHROM_SIZES} > non_exons.bed
    bedtools complement -i ${INFOLDER}/parsed_mRNA_intron.merged.bed -g ${CHROM_SIZES} > non_introns.bed

    ## promoters > exon
    ## first copy the input promoter file into the final output folder, to be consistent
    cp ${INFOLDER}/parsed_mRNA_promoters.merged.bed final_files/parsed_mRNA_promoters.merged.bed
    bedtools intersect -a ${INFOLDER}/parsed_mRNA_exon.merged.bed -b non_promoters.bed > final_files/np_exons.bed

    ## exon > intron
    bedtools intersect -a ${INFOLDER}/parsed_mRNA_intron.merged.bed -b non_promoters.bed > np_introns.bed
    bedtools intersect -a np_introns.bed -b non_exons.bed > final_files/npe_introns.bed

    ## intron > repeats
    bedtools intersect -a ${INFOLDER}/parsed_repeats.merged.bed -b non_promoters.bed > np_repeats.bed
    bedtools intersect -a np_repeats.bed -b non_exons.bed > npe_repeats.bed
    bedtools intersect -a npe_repeats.bed -b non_introns.bed > final_files/npei_repeats.bed

    echo "Sorting output files"

    for F in final_files/*.bed; do
	echo $F
	sort -k1,1V -k2,2n $F > $F.sorted
	mv $F.sorted $F
    done
    
    cd -
    
else
    echo "Usage: $0 <chromosome sizes file> <input folder> <output folder>"
fi
