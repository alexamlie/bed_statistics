"""
generate_isoforms.py
Alex Amlie-Wolf, Spring 2014

A program that takes in a refseq formatted file containing exon starts and ends for different
transcripts and outputs a bed file containing windows for all of the exons and introns of that
gene. This can't be done (as far as I know) directly from the UCSC table downloads or bedtools,
as they restrict you to making a bed interval for either exons or introns, or the whole gene, etc

Only works with Python 2.7 and above
"""

import argparse

def ref_to_bed(ref_file, out_file):
    # define a convenience dict to access the refseq file
    ref_coords = {'bin': 0, 'cdsEnd': 7, 'cdsEndStat': 14, 'cdsStart': 6, 'cdsStartStat': 13,
 'chrom': 2, 'exonCount': 8, 'exonEnds': 10, 'exonFrames': 15, 'exonStarts': 9, 'name': 1,
 'name2': 12, 'score': 11, 'strand': 3, 'txEnd': 5, 'txStart': 4}

    with open(ref_file, 'r') as refseq, open(out_file, 'w') as outfile:
        for gene in refseq:
            genedata = gene.strip().split("\t")

            # to go through the windows of the gene (introns and exons), we need to make them
            # into lists. skips the last character because it's a comma
            exonstarts = genedata[ref_coords['exonStarts']][:-1].split(',')
            exonends = genedata[ref_coords['exonEnds']][:-1].split(',')

            # build up the lists of introns
            intronstarts = []
            intronends = []
            # first see if the gene starts with an intron (may not ever happen)
            if int(exonstarts[0]) > int(genedata[ref_coords["txStart"]]):
                intronstarts.append(genedata[ref_coords["txStart"]])
                intronends.append(exonstarts[0])
            # loop through exons (this assumes the lists are the same length)
            for i in range(len(exonstarts)-1):
                intronstarts.append(exonends[i]) # intron starts where exon ends
                intronends.append(exonstarts[i+1]) # intron ends where exon starts
            # now check if there is a 3' intron/UTR/whatever
            if int(exonends[len(exonends)-1]) < int(genedata[ref_coords["txEnd"]]):
                intronstarts.append(exonends[len(exonends)-1])
                intronends.append(genedata[ref_coords['txEnd']])

            # now just make the bed lines
            # bed format: chr start end name score strand
            ## write the exons
            for i in range(len(exonstarts)):
                outfile.write("\t".join([genedata[ref_coords['chrom']], exonstarts[i], exonends[i], genedata[ref_coords["name"]]+"_exon_"+str(i), "0", genedata[ref_coords["strand"]]])+"\n")

            ## write the introns
            for i in range(len(intronstarts)):
                outfile.write("\t".join([genedata[ref_coords['chrom']], intronstarts[i], intronends[i], genedata[ref_coords["name"]]+"_intron_"+str(i), "0", genedata[ref_coords["strand"]]])+"\n")

        
if __name__=="__main__":
    # create argparser
    parser = argparse.ArgumentParser(description="Turn a refseq file into a .bed file with exon and intron windows")
    parser.add_argument("refseq_file", help="The refseq file. Should be in refseq format")
    parser.add_argument("output", help="The output file you wish to use")
    pargs = parser.parse_args()
    ref_to_bed(pargs.refseq_file, pargs.output)
