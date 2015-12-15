#!/bin/bash

## generate_unstranded_elements.sh
## alex amlie-wolf
## instead of using stranded partitions, this script looks only at the collapsed elements
## optional fourth element tells it to use UTRs and describes how you want it

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
    
elif [ $# == 4 ]; then
    UTR_message=$1
    if [ "${UTR_message}" == "UTR_exons_introns" ]; then
	echo "Parsing unstranded partition including UTR exons and introns"
	CHROM_SIZES=$2
	INFOLDER=$3
	OUTFOLDER=$4
	
	## assumes that the input folder is organized the same way as the
	## 12_04_15_protein_coding_file_generation.sh script
	mkdir -p ${OUTFOLDER}/final_files/ 
	
	cd ${OUTFOLDER}/

	## generate all the complements:	
	bedtools complement -i ${INFOLDER}/parsed_UTR5_exon.merged.bed -g ${CHROM_SIZES} > non_5utr_exons.bed
	bedtools complement -i ${INFOLDER}/parsed_UTR5_intron.merged.bed -g ${CHROM_SIZES} > non_5utr_introns.bed
	bedtools complement -i ${INFOLDER}/parsed_UTR3_exon.merged.bed -g ${CHROM_SIZES} > non_3utr_exons.bed
	bedtools complement -i ${INFOLDER}/parsed_UTR3_intron.merged.bed -g ${CHROM_SIZES} > non_3utr_introns.bed
	bedtools complement -i ${INFOLDER}/parsed_mRNA_promoters.merged.bed -g ${CHROM_SIZES} > non_promoters.bed
	bedtools complement -i ${INFOLDER}/parsed_mRNA_exon.merged.bed -g ${CHROM_SIZES} > non_exons.bed
	bedtools complement -i ${INFOLDER}/parsed_mRNA_intron.merged.bed -g ${CHROM_SIZES} > non_introns.bed

	## 5' UTR exon > 5' UTR intron
	## first copy the input 5' UTR exon file into the final output folder, to be consistent
	cp ${INFOLDER}/parsed_UTR5_exon.merged.bed final_files/parsed_5utr_exons.merged.bed	
	bedtools intersect -a ${INFOLDER}/parsed_UTR5_intron.merged.bed -b non_5utr_exons.bed > final_files/n5e_5utr_introns.bed

	## 5' UTR intron > 3' UTR exon
	bedtools intersect -a ${INFOLDER}/parsed_UTR3_exon.merged.bed -b non_5utr_exons.bed > n5e_3utr_exons.bed
	bedtools intersect -a n5e_3utr_exons.bed -b non_5utr_introns.bed > final_files/n5e5i_3utr_exons.bed

	## 3' UTR exon > 3' UTR intron
	bedtools intersect -a ${INFOLDER}/parsed_UTR3_intron.merged.bed -b non_5utr_exons.bed > n5e_3utr_introns.bed
	bedtools intersect -a n5e_3utr_introns.bed -b non_5utr_introns.bed > n5e5i_3utr_introns.bed
	bedtools intersect -a n5e5i_3utr_introns.bed -b non_3utr_exons.bed > final_files/n5e5i3e_3utr_introns.bed

	## 3' UTR intron > promoter
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_promoters.merged.bed -b non_5utr_exons.bed > n5e_promoters.bed
	bedtools intersect -a n5e_promoters.bed -b non_5utr_introns.bed > n5e5i_promoters.bed
	bedtools intersect -a n5e5i_promoters.bed -b non_3utr_exons.bed > n5e5i3e_promoters.bed
	bedtools intersect -a n5e5i3e_promoters.bed -b non_3utr_introns.bed > final_files/n5e5i3e3i_promoters.bed

	## promoter > exon
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_exon.merged.bed -b non_5utr_exons.bed > n5e_exons.bed
	bedtools intersect -a n5e_exons.bed -b non_5utr_introns.bed > n5e5i_exons.bed
	bedtools intersect -a n5e5i_exons.bed -b non_3utr_exons.bed > n5e5i3e_exons.bed
	bedtools intersect -a n5e5i3e_exons.bed -b non_3utr_introns.bed > n5e5i3e3i_exons.bed
	bedtools intersect -a n5e5i3e3i_exons.bed -b non_promoters.bed > final_files/n5e5i3e3ip_exons.bed

	## exon > intron
	bedtools intersect -a ${INFOLDER}/parsed_mRNA_intron.merged.bed -b non_5utr_exons.bed > n5e_introns.bed
	bedtools intersect -a n5e_introns.bed -b non_5utr_introns.bed > n5e5i_introns.bed
	bedtools intersect -a n5e5i_introns.bed -b non_3utr_exons.bed > n5e5i3e_introns.bed
	bedtools intersect -a n5e5i3e_introns.bed -b non_3utr_introns.bed > n5e5i3e3i_introns.bed
	bedtools intersect -a n5e5i3e3i_introns.bed -b non_promoters.bed > n5e5i3e3ip_introns.bed
	bedtools intersect -a n5e5i3e3ip_introns.bed -b non_exons.bed > final_files/n5e5i3e3ipe_introns.bed

	## intron > repeat
	bedtools intersect -a ${INFOLDER}/parsed_repeats.merged.bed -b non_5utr_exons.bed > n5e_repeats.bed
	bedtools intersect -a n5e_repeats.bed -b non_5utr_introns.bed > n5e5i_repeats.bed
	bedtools intersect -a n5e5i_repeats.bed -b non_3utr_exons.bed > n5e5i3e_repeats.bed
	bedtools intersect -a n5e5i3e_repeats.bed -b non_3utr_introns.bed > n5e5i3e3i_repeats.bed
	bedtools intersect -a n5e5i3e3i_repeats.bed -b non_promoters.bed > n5e5i3e3ip_repeats.bed
	bedtools intersect -a n5e5i3e3ip_repeats.bed -b non_exons.bed > n5e5i3e3ipe_repeats.bed
	bedtools intersect -a n5e5i3e3ipe_repeats.bed -b non_introns.bed > final_files/n5e5i3e3ipei_repeats.bed

	echo "Sorting output files"

	for F in final_files/*.bed; do
	    echo $F
	    sort -k1,1V -k2,2n $F > $F.sorted
	    mv $F.sorted $F
	done
	
	cd -

    elif [ "${UTR_message}" == "UTRs" ]; then
	echo "Non-exon/intron analysis not currently supported"
    else
	echo "UTR argument must equal either 'UTR_exons_introns' or 'UTRs'"
    fi
else
    echo "Usage: $0 [UTR_exons_introns | UTRs ] <chromosome sizes file> <input folder> <output folder>"
    echo "First argument is optional and tells the script to use UTR annotations, choosing between using introns and exons, or just using full UTR annotations"
fi
