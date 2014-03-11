"""
bed_statistics.py
Alex Amlie-Wolf, Spring 2014

A program to compute various coverage statistics given a .bed file and a reference genome
Only works with Python 2.7 and above, because of the with open, open, open construction
"""

import argparse, sys, os

def compute_statistics(reffile, bedfile, output):
    # define convenience dicts to access the refseq file and the bedfile
    ref_coords = {'bin': 0, 'cdsEnd': 7, 'cdsEndStat': 14, 'cdsStart': 6, 'cdsStartStat': 13,
 'chrom': 2, 'exonCount': 8, 'exonEnds': 10, 'exonFrames': 15, 'exonStarts': 9, 'name': 1,
 'name2': 12, 'score': 11, 'strand': 3, 'txEnd': 5, 'txStart': 4}

    bed_coords = {} # this is standard for .bed files
    bed_coords["chrom"] = 0
    bed_coords["start"] = 1
    bed_coords["end"] = 2
    bed_coords["strand"] = 5
    
    with open(reffile, 'r') as refseq, open(bedfile, 'r') as bedpeaks, open(output, 'w') as outfile:
        # list of statistcs to compute
        # refseq name, normal name, number of exons hit, total number of exons, proportion of
        # exonic bases hit, number of exonic bases covered, same 4 for introns, total bases
        # covered, proportion of total bases covered, and total number of bases
        stats = ('refseq', 'symbol', 'coveredExons', 'totalExons', 'propExonBases',
        'exonBasesCovered', 'exonicHits', 'coveredIntrons', 'totalIntrons', 'propIntronBases',
        'intronBasesCovered', 'intronicHits', 'totalBasesCovered', 'propBasesCovered',
        'totalBases')
        outfile.write("\t".join(stats)+"\n") # write the header

        # we don't have to parse any of the refseq or bed file lines because we assume
        # that pre-processing has been done to sort them and remove the headers
        
        # this data structure will hold all the peaks we are currently considering
        # it gets updated as we go through the genes to get rid of peaks that can no
        # longer possibly match. it stores them as arrays for ease
        curpks = []

        # initialize the peak tracking
        curpk = bedpeaks.readline().strip()
        curpkdata = curpk.split("\t")
        curpks.append(curpkdata)
        
        # now iterate through the refseq genes
        for gene in refseq:
            # make this gene into an array so we can use our ref_coords dict
            genedata = gene.strip().split("\t")
            # initialize the output string for this gene
            outstring = "\t".join([genedata[ref_coords["name"]], genedata[ref_coords["name2"]]])

            # we now must get rid of each peak that cannot overlap with this gene
            for i in range(len(curpks)-1, -1, -1): # loop backwards to remove
                pk = curpks[i]
                # if we end before the gene starts, are on the wrong chromosome, or the wrong
                # strand, we no longer consider the peak
                # ( int(pk[bed_coords['end']]) < int(genedata[ref_coords['txStart']]) ) or
                if (pk[bed_coords['chrom']] != genedata[ref_coords['chrom']]) or (pk[bed_coords['strand']] != genedata[ref_coords['strand']]):
                    curpks.pop(i) # remove this peak

            # now, we need to add all the possible next peaks that might match
            # we add peaks until we violate one of the conditions seen above
            # also, make sure that we haven't finished reading through the file
            # lengthbool = int(curpkdata[bed_coords['end']]) >= int(genedata[ref_coords['txStart']]) and int(curpkdata[bed_coords['start']]) <= int(genedata[ref_coords['txEnd']])
            chrombool = curpkdata[bed_coords['chrom']] == genedata[ref_coords['chrom']]
            strandbool = curpkdata[bed_coords['strand']] == genedata[ref_coords['strand']] 
            # while curpk and lengthbool and chrombool and strandbool:
            while curpk and chrombool and strandbool:
                curpks.append(curpkdata) # add the current peak
                curpk = bedpeaks.readline().strip() # get the next peak to check
                curpkdata = curpk.split("\t") # split it to get its juicy data
                # lengthbool = int(curpkdata[bed_coords['end']]) >= int(genedata[ref_coords['txStart']]) and int(curpkdata[bed_coords['start']]) <= int(genedata[ref_coords['txEnd']])
                chrombool = curpkdata[bed_coords['chrom']] == genedata[ref_coords['chrom']]
                strandbool = curpkdata[bed_coords['strand']] == genedata[ref_coords['strand']] 

            print "Curpk number currently", len(curpks)
                                
            # now we have all the possible overlapping peaks in the curpks list
            exons_covered = 0 # track the number of exons covered
            exonbases = 0 # track the number of exonic bases covered
            exonhits = 0 # track the number of peaks overlapping exons
            introns_covered = 0 # track the number of introns covered
            intronbases = 0 # track the number of intronic bases covered
            intronhits = 0 # track the number of peaks overlapping introns

            # to go through the windows of the gene (introns and exons), we need to make them
            # into lists. skips the last character because it's a comma
            exonstarts = map(int, genedata[ref_coords['exonStarts']][:-1].split(','))
            exonends = map(int, genedata[ref_coords['exonEnds']][:-1].split(','))

            # build up the lists of introns
            intronstarts = []
            intronends = []
            # first see if the gene starts with an intron (may not ever happen)
            if exonstarts[0] > int(genedata[ref_coords["txStart"]]):
                intronstarts.append(int(genedata[ref_coords["txStart"]]))
                intronends.append(exonstarts[0])
            # loop through exons (this assumes the lists are the same length)
            for i in range(len(exonstarts)-1):
                intronstarts.append(exonends[i]) # intron starts where exon ends
                intronends.append(exonstarts[i+1]) # intron ends where exon starts
            # now check if there is a 3' intron/UTR/whatever
            if exonends[len(exonends)-1] < int(genedata[ref_coords["txEnd"]]):
                intronstarts.append(exonends[len(exonends)-1])
                intronends.append(int(genedata[ref_coords['txEnd']]))

            # now we have our information on where the exons and introns are, so we do analysis
            # loop through exons
            for i in range(len(exonstarts)):
                covered = False
                start = exonstarts[i]
                end = exonends[i]
                intervalsCovered = [] # this tracks which bases have been covered
                for pk in curpks:
                    # first figure out if we've gone too far past this exon, to avoid
                    # unecessary looping
                    if int(pk[bed_coords['start']]) > end:
                        continue # try checking more peaks?
                    # then check if this peak can even overlap (we have to do this to try
                    # to optimize because we only remove peaks from curpks at each gene)
                    if int(pk[bed_coords['end']]) < start:
                        continue # skip this peak
                    # if we reach this point, then our peak covers at least part of this exon
                    covered = True
                    exonhits += 1 # add another exonic peak
                    
                    # now compute the coverage
                    interval = [max(start, int(pk[bed_coords["start"]])), min(end, int(pk[bed_coords["end"]]))]
                    # check if this interval overlaps with any in the list we have so far
                    overlap = False
                    for i in range(len(intervalsCovered)):
                        curint = intervalsCovered[i]
                        # check if our end is after their start and our start is before
                        # their end. this will catch both partial and full overlaps
                        if interval[1] >= curint[0] and interval[0] <= curint[1]:
                            overlap = True
                            # expand the interval
                            intervalsCovered[i] = [min(interval[0], curint[0]), max(interval[1], curint[1])]

                    # if we didn't overlap, add it to the list
                    if not overlap:
                        intervalsCovered.append(interval)
                        
                if covered:
                    exons_covered += 1
                    # we first need to process our intervals list to merge overlaps
                    mergedIntervals = mergeOverlaps(intervalsCovered)
                    # now we get all the bases covered
                    for interval in mergedIntervals:
                        exonbases += interval[1] - interval[0]

            # done looping through exons
            # calculate the proportion of exonic bases covered
            exon_prop = (float(exonbases) / float(sum([end-start for (end, start) in zip(exonends, exonstarts)]))) if exonbases > 0 else 0
            # add the info to the output string: number of exons covered, total exons
            # proportion of exonic bases, number of exonic bases covered, and number of exonic
            # peaks
            outstring = "\t".join([outstring, str(exons_covered), str(len(exonstarts)), str(exon_prop), str(exonbases), str(exonhits)])

            # loop through introns now
            for i in range(len(intronstarts)):
                covered = False
                start = intronstarts[i]
                end = intronends[i]
                intervalsCovered = [] # this tracks which bases have been covered
                for pk in curpks:
                    # first figure out if we've gone too far past this intron, to avoid
                    # unecessary looping
                    if int(pk[bed_coords['start']]) > end:
                        break
                    # then check if this peak can even overlap (we have to do this to try
                    # to optimize because we only remove peaks from curpks at each gene)
                    if int(pk[bed_coords['end']]) < start:
                        continue # skip this peak
                    # if we reach this point, then our peak covers at least part of this intron
                    covered = True
                    intronhits += 1 # add another intronic peak
                    
                    # now compute the coverage
                    interval = [max(start, int(pk[bed_coords["start"]])), min(end, int(pk[bed_coords["end"]]))]
                    # check if this interval overlaps with any in the list we have so far
                    overlap = False
                    for i in range(len(intervalsCovered)):
                        curint = intervalsCovered[i]
                        # check if our end is after their start and our start is before
                        # their end. this will catch both partial and full overlaps
                        # note that this just expands each interval individually
                        # the merging step is done later
                        if interval[1] >= curint[0] and interval[0] <= curint[1]:
                            overlap = True
                            # expand the interval
                            intervalsCovered[i] = [min(interval[0], curint[0]), max(interval[1], curint[1])]

                    # if we didn't overlap, add it to the list
                    if not overlap:
                        intervalsCovered.append(interval)
                        
                if covered:
                    introns_covered += 1
                    # we first need to process our intervals list to merge overlaps
                    mergedIntervals = mergeOverlaps(intervalsCovered)
                    # now we get all the bases covered
                    for interval in mergedIntervals:
                        intronbases += interval[1] - interval[0]

            # done looping through introns
            # calculate the proportion of intronic bases covered
            intron_prop = (float(intronbases) / float(sum([end-start for (end, start) in zip(intronends, intronstarts)]))) if intronbases > 0 else 0
            # add the info to the output string: number of introns covered, total introns
            # proportion of intronic bases, number of intronic bases covered, # of intronic peaks
            outstring = "\t".join([outstring, str(introns_covered), str(len(intronstarts)), str(intron_prop), str(intronbases), str(intronhits)])

            # finally, add the remaining statistics to the output string and write it
            # calculate the total proportion of bases covered
            total_prop = (float(intronbases+exonbases) / (float(genedata[ref_coords['txEnd']]) - float(genedata[ref_coords['txStart']]))) if intronbases+exonbases > 0 else 0
            outstring = "\t".join([outstring, str(intronbases+exonbases), str(total_prop), str(float(genedata[ref_coords['txEnd']]) - float(genedata[ref_coords['txStart']]))])
            outfile.write(outstring+"\n")

        print "Analysis complete!"
                 
def mergeOverlaps(intervalsCovered):
    retints = [intervalsCovered[0]] # start with the first interval
    for start, end in intervalsCovered:
        if start < retints[-1][1]: # if we start before the most recent interval ends
            retints[-1][1] = max(end, retints[-1][1]) # expand the interval
        else: # no overlap
            retints.append([start, end])
    return retints

if __name__=="__main__":
    # create the argument parser
    parser = argparse.ArgumentParser(description="Compute coverage statistics of a .bed file over refseq")
    parser.add_argument("reffile", help="The sorted reference file. Should be in the refseq format, and sorted according to chromosome, then transcript start.")
    parser.add_argument("bedfile", help="The .bed file that you want to compute overlap statistics for.")
    parser.add_argument("output", help="The name of the file you want to write the statistics to.")
    pargs = parser.parse_args()
    compute_statistics(pargs.reffile, pargs.bedfile, pargs.output)

