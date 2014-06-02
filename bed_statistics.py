"""
bed_statistics.py
Alex Amlie-Wolf, Spring 2014

A program to compute various coverage statistics given a .bed file and a reference
genome in refseq format.
Only works with Python 2.7 and above, because of the with open, open, open construction
Both the .bed and the reference file must be sorted by chromosome and strand, then by start
position. This is to ensure that we don't skip any peaks that we should look at
"""

import argparse, sys, os

"""
Takes in the reference file, the dict telling us how to refer to that reference file,
the input .bed file, the output file, and the log file
"""
def compute_statistics(reffile, ref_dict, bedfile, output, log):
    # define convenience dicts to access the refseq file and the bedfile
    bed_coords = {} 
    bed_coords["chrom"] = 0
    bed_coords["start"] = 1
    bed_coords["end"] = 2
    bed_coords["strand"] = 5

    # standardize the reference dict (refseq and common name, transcription start, etc)
    # we want refseq to be 'name' and the gene symbol to be 'name2'
    if 'refseq' in ref_dict.keys():
        ref_dict['name'] = ref_dict.pop('refseq')
    if 'geneSymbol' in ref_dict.keys():
        ref_dict['name2'] = ref_dict.pop('geneSymbol')

    # we want the start of the transcript to be txStart and the end to be txEnd
    if 'chromStart' in ref_dict.keys():
        ref_dict['txStart'] = ref_dict.pop('chromStart')
    if 'chromEnd' in ref_dict.keys():
        ref_dict['txEnd'] = ref_dict.pop('chromEnd')    
        
        
    with open(reffile, 'r') as reference, open(bedfile, 'r') as bedpeaks, open(output, 'w') as outfile, open(log, 'w') as logfile:
        # list of statistcs to compute refseq name, normal name, number of exons hit, total
        # number of exons, proportion of exonic bases hit, number of exonic bases covered, same
        # 4 for introns, same 4 for 3' and 5' UTRs, total bases covered, proportion of total
        # bases covered, total number of bases, and total number of CLIP hits (for convenience)
        stats = ('refseq', 'symbol', 'covered_exons', 'total_exons', 'prop_exon_bases',
        'exon_bases_covered', 'exonic_bases', 'exonic_hits', 'covered_introns',
        'total_introns', 'prop_intron_bases', 'intron_bases_covered', 'intronic_bases',
        'intronic_hits', 'tp_utr_bases', 'tp_utr_bases_covered', 'tp_utr_hits', 'tp_utr_prop',
        'fp_utr_bases', 'fp_utr_bases_covered', 'fp_utr_hits', 'fp_utr_prop',
        'total_bases_covered', 'prop_bases_covered', 'total_bases', 'total_hits')
        outfile.write("\t".join(stats)+"\n") # write the header

        # we don't have to pre-parse any of the reference or bed file lines because we assume
        # that pre-processing has been done to sort them and remove the headers
        
        # this data structure will hold all the peaks we are currently considering
        # it gets updated as we go through the genes to get rid of peaks that can no
        # longer possibly match. it stores them as arrays for ease
        curpks = []

        # initialize the peak tracking variables
        # what i'm doing here is making sure that the first iteration of the loop
        # goes through. These all get reset as soon as the script starts
        curpk = bedpeaks.readline().strip()
        curpkdata = curpk.split("\t") # don't append this data yet
        
        # now iterate through the reference genes
        for gene in reference:            
            # make this gene into an array so we can use our ref_coords dict
            genedata = gene.strip().split("\t")
            # initialize the output string for this gene            
            outstring = "\t".join([genedata[ref_dict["name"]], genedata[ref_dict["name2"]]])
            logfile.write("Computing statistics on "+"\t".join([genedata[ref_dict["name"]], genedata[ref_dict["name2"]], genedata[ref_dict["chrom"]], genedata[ref_dict["txStart"]], genedata[ref_dict["txEnd"]], genedata[ref_dict['strand']]])+"\n")
                        
            # we now must get rid of each peak that cannot overlap with this gene
            # this is to reduce memory use
            for i in range(len(curpks)-1, -1, -1): # loop backwards to remove
                pk = curpks[i]
                # if we are on the wrong chromosome or the wrong strand, we no longer consider
                # the peak
                # ( int(pk[bed_coords['end']]) < int(genedata[ref_dict['txStart']]) ) or
                if (pk[bed_coords['chrom']] != genedata[ref_dict['chrom']]) or (pk[bed_coords['strand']] != genedata[ref_dict['strand']]):
                    curpks.pop(i) # remove this peak

            # now, we need to add all the possible next peaks that might match
            # we add peaks until we violate one of the conditions seen above
            # also, make sure that we haven't finished reading through the file
            # lengthbool = int(curpkdata[bed_coords['end']]) >= int(genedata[ref_dict['txStart']]) and int(curpkdata[bed_coords['start']]) <= int(genedata[ref_dict['txEnd']])
            # here, we check curpk first so that if it's false, we don't try to evaluate
            # curpkdata and get an index out of bounds error
            chrombool = curpk and curpkdata[bed_coords['chrom']] == genedata[ref_dict['chrom']]
            strandbool = curpk and curpkdata[bed_coords['strand']] == genedata[ref_dict['strand']]
       
            while chrombool and strandbool:
                curpks.append(curpkdata) # append the current peak
                curpk = bedpeaks.readline().strip() # get the next peak to check
                if curpk: # only if we're not done reading the file
                    curpkdata = curpk.split("\t") # split it to get its juicy data
                    # lengthbool = int(curpkdata[bed_coords['end']]) >= int(genedata[ref_dict['txStart']]) and int(curpkdata[bed_coords['start']]) <= int(genedata[ref_dict['txEnd']])
                    chrombool = curpkdata[bed_coords['chrom']] == genedata[ref_dict['chrom']]
                    strandbool = curpkdata[bed_coords['strand']] == genedata[ref_dict['strand']]
                else: 
                    break
                                                    
            # now we have all the possible overlapping peaks in the curpks list

            # we have 4 categories of gene section:
            # 5' UTR (first exon start, which is the same as txStart -> cdsStart)
            # coding exons (cdsStart->first exon end, second exon start->second exon end, ..., last exon start -> cdsEnd)
            # introns (second exon start-first exon end, third exon start-second exon end, ...)
            # 3' UTR (cdsEnd -> last exon end, which is the same as txEnd)
            
            exons_covered = 0 # track the number of exons covered
            exonbases = 0 # track the number of exonic bases covered
            exonhits = 0 # track the number of peaks overlapping exons
            introns_covered = 0 # track the number of introns covered
            intronbases = 0 # track the number of intronic bases covered
            intronhits = 0 # track the number of peaks overlapping introns
            tp_utrbases = 0 # track the 3' utr bases covered
            tp_utrhits = 0 # track the number of 3' UTR hits
            fp_utrbases = 0 # track the 5' UTR bases covered
            fp_utrhits = 0 # track the number of 5' UTR hits

            # define the 4 categories
            # start with the UTRs, which depend on strand
            if genedata[ref_dict['strand']] == '+':
                # in the plus direction, 5' starts at the beginning of the
                # transcript, 3' ends at the end
                fp_utr_start = int(genedata[ref_dict['txStart']])
                fp_utr_end = int(genedata[ref_dict['cdsStart']])
                tp_utr_start = int(genedata[ref_dict['cdsEnd']])
                tp_utr_end = int(genedata[ref_dict['txEnd']])
            else:
                # in the minus direction, 5' starts at the end, 3' ends
                # at the beginning
                # this is a bit confusing, because the 'starts' of the UTRs as defined by
                # these variables go in the 'plus' direction. this is so that we can make the
                # same index-based comparisons regardless of whether we're on the + or - strand
                # i.e. the 5' UTR will actually go from the end of the transcript to the end
                # of the coding region on the minus strand ( cdsEnd <-------- txEnd ), but
                # we define it the other way so we still ask if a peak ends before it starts
                # (functionally, cdsEnd ---------> txEnd )
                fp_utr_start = int(genedata[ref_dict['cdsEnd']])
                fp_utr_end = int(genedata[ref_dict['txEnd']])
                tp_utr_start = int(genedata[ref_dict['txStart']])
                tp_utr_end = int(genedata[ref_dict['cdsStart']])

            # now define the coding exons
            # if the coding region start and end are the same, we have no coding exons
            if int(genedata[ref_dict['cdsStart']]) == int(genedata[ref_dict['cdsEnd']]):
                logfile.write("Non-coding transcript\n")
                outstring = "\t".join([outstring, "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"])
            else:
                # start by using all the exons
                # skips the last character because it's a comma
                exonstarts = map(int, genedata[ref_dict['exonStarts']][:-1].split(','))
                exonends = map(int, genedata[ref_dict['exonEnds']][:-1].split(','))
                # however, these lists contain the UTRs, because of the way refseq works
                # (the first exon starts at the beginning of the transcript and the last exon
                # ends at the end of the transcript, always), so we have to change them

                # first coding exon starts at the beginning of the coding region            
                # remove exons that are in the wrong place
                # loop backwards so that we can remove
                for i in range(len(exonends)-1, -1, -1):
                    # if an exon ends before the coding region begins, or starts after the
                    # coding region ends, remove it and go to the next one
                    if exonends[i] < int(genedata[ref_dict['cdsStart']]) or exonstarts[i] > int(genedata[ref_dict['cdsEnd']]):
                        exonends.pop(i)
                        exonstarts.pop(i)
                        continue

                    # check for whether we need to set this start to the cdsStart
                    # this is when the exon starts before the cdsEnd and ends after cdsStart
                    # (which we already know, from the test above)
                    if exonstarts[i] < int(genedata[ref_dict['cdsStart']]):
                        exonstarts[i] = int(genedata[ref_dict['cdsStart']])
                    # check for whether the end needs to be set to cdsEnd
                    # this is when the start is before cdsEnd (which we already know it is),
                    # and the end is after it
                    if exonends[i] > int(genedata[ref_dict['cdsEnd']]):
                        exonends[i] = int(genedata[ref_dict['cdsEnd']])
                
                # build up the lists of introns
                intronstarts = []
                intronends = []
                # loop through exons (this assumes the two exon lists are the same length,
                # which they always should be)
                for i in range(len(exonstarts)-1):
                    intronstarts.append(exonends[i]) # intron starts where exon ends
                    intronends.append(exonstarts[i+1]) # intron ends where exon starts

                # now we have our information on where the exons and introns are, so we do
                # analysis print out the lists to the log file
                logfile.write("\t".join(["5' UTR:", str(fp_utr_start), str(fp_utr_end), "exons:", str([[exonst, exonend] for (exonst, exonend) in zip(exonstarts, exonends)]), "introns:", str([[intronst, intronend] for (intronst, intronend) in zip(intronstarts, intronends)]), "3' UTR:", str(tp_utr_start), str(tp_utr_end)])+"\n")

                # -------------------------------------------------
                # loop through exons
                # -------------------------------------------------
                for i in range(len(exonstarts)):
                    covered = False
                    start = exonstarts[i]
                    end = exonends[i]
                    intervalsCovered = [] # this tracks which bases have been covered
                    for pk in curpks:                                        
                        # first figure out if we've gone too far past this exon, to avoid
                        # unecessary looping
                        if int(pk[bed_coords['start']]) > end:
                            break # we're too far, so stop checking
                        # then check if this peak can possibly overlap 
                        # (we only remove peaks from curpks at each gene iteration)
                        if int(pk[bed_coords['end']]) < start:
                            continue # skip this peak
                        # if we reach this point, then our peak covers at least part of this
                        # exon
                        logfile.write("Peak overlap: peak "+"\t".join([pk[bed_coords['start']], pk[bed_coords['end']]])+" exon "+"\t".join([str(start), str(end)])+"\n")
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
                        # mergedIntervals = mergeOverlaps(intervalsCovered)
                        # now we get all the bases covered
                        for interval in intervalsCovered:
                            exonbases += interval[1] - interval[0]

                # done looping through exons
                # calculate the number of exonic bases
                exon_base_total = float(sum([end-start for (end, start) in zip(exonends, exonstarts)]))
                # calculate the proportion of exonic bases covered
                exon_prop = (float(exonbases) / exon_base_total) if exonbases > 0 else 0
                # add the info to the output string: number of exons covered, total exons
                # proportion of exonic bases, number of exonic bases covered, and number of exonic
                # peaks
                outstring = "\t".join([outstring, str(exons_covered), str(len(exonstarts)), str(exon_prop), str(exonbases), str(exon_base_total), str(exonhits)])

                # -------------------------------------------------            
                # loop through introns now
                # -------------------------------------------------            
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
                        ## mergedIntervals = mergeOverlaps(intervalsCovered)
                        # now we get all the bases covered
                        for interval in intervalsCovered:
                            intronbases += interval[1] - interval[0]

                # done looping through introns
                # calculate the number of intronic bases
                intron_base_total = float(sum([end-start for (end, start) in zip(intronends, intronstarts)]))

                # calculate the proportion of intronic bases covered
                intron_prop = (float(intronbases) / intron_base_total) if intronbases > 0 else 0
                # add the info to the output string: number of introns covered, total introns
                # proportion of intronic bases, number of intronic bases covered, # of intronic peaks
                outstring = "\t".join([outstring, str(introns_covered), str(len(intronstarts)), str(intron_prop), str(intronbases), str(intron_base_total), str(intronhits)])

            # this ends the indentation block for whether there is a coding region or not
            # -------------------------------------------------            
            # now do 3' UTRs
            # -------------------------------------------------            
            intervalsCovered = [] 
            covered = False
            for pk in curpks:
                # if this peak starts after the UTR ends, we stop looking
                if int(pk[bed_coords['start']]) > tp_utr_end:
                    break
                # if it ends before this UTR begins, we skip it
                if int(pk[bed_coords['end']]) < tp_utr_start:
                    continue
                # now our peak overlaps
                logfile.write("Peak overlap of 3' UTR: "+"\t".join([pk[bed_coords['start']], pk[bed_coords['end']]])+" UTR region "+"\t".join([str(tp_utr_start), str(tp_utr_end)])+"\n")
                covered = True
                tp_utrhits += 1 # add a 3' UTR hit
                interval = [max(tp_utr_start, int(pk[bed_coords['start']])), min(tp_utr_end, int(pk[bed_coords['end']]))]
                
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

            # done looping through, compute the bases covered
            if covered:
                # mergedIntervals = mergeOverlaps(intervalsCovered)
                for interval in intervalsCovered:
                    tp_utrbases += interval[1] - interval[0]

            tp_utr_base_total = tp_utr_end - tp_utr_start
            
            # add the 3' UTR information
            tp_prop = (float(tp_utrbases)/float(tp_utr_base_total)) if tp_utr_base_total > 0 else 0
            outstring = "\t".join([outstring, str(tp_utr_base_total), str(tp_utrbases), str(tp_utrhits), str(tp_prop)])

            # -------------------------------------------------                        
            # now do 5' UTRs
            # -------------------------------------------------            
            intervalsCovered = [] 
            covered = False
            for pk in curpks:
                # if this peak starts after the UTR ends, we stop looking
                if int(pk[bed_coords['start']]) > fp_utr_end:
                    break
                # if it ends before this UTR begins, we skip it
                if int(pk[bed_coords['end']]) < fp_utr_start:
                    continue
                # now our peak overlaps
                logfile.write("Peak overlap of 5' UTR: "+"\t".join([pk[bed_coords['start']], pk[bed_coords['end']]])+" UTR region "+"\t".join([str(fp_utr_start), str(fp_utr_end)])+"\n")
                covered = True
                fp_utrhits += 1 # add a 5' UTR hit
                interval = [max(fp_utr_start, int(pk[bed_coords['start']])), min(fp_utr_end, int(pk[bed_coords['end']]))]
                
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

            # done looping through, compute the bases covered
            if covered:
                ## mergedIntervals = mergeOverlaps(intervalsCovered)
                for interval in intervalsCovered:
                    fp_utrbases += interval[1] - interval[0]

            fp_utr_base_total = fp_utr_end - fp_utr_start
            
            # add the 5' UTR information
            fp_prop = (float(fp_utrbases)/float(fp_utr_base_total)) if fp_utr_base_total > 0 else 0
            outstring = "\t".join([outstring, str(fp_utr_base_total), str(fp_utrbases), str(fp_utrhits), str(fp_prop)])
                        
            # finally, add the remaining statistics to the output string and write it
            # calculate the total proportion of bases covered
            total_prop = (float(intronbases+exonbases+tp_utrbases+fp_utrbases) / (float(genedata[ref_dict['txEnd']]) - float(genedata[ref_dict['txStart']]))) if intronbases+exonbases+tp_utrbases_fp_utrbases > 0 else 0
            outstring = "\t".join([outstring, str(intronbases+exonbases+tp_utrbases+fp_utrbases), str(total_prop), str(float(genedata[ref_dict['txEnd']]) - float(genedata[ref_dict['txStart']])), str(exonhits+intronhits+tp_utrhits+fp_utrhits)])
            
            outfile.write(outstring+"\n")

        print "Analysis complete!"

'''
This method takes a list of intervals and merges any overlapping ones together.
I'm leaving it in for completeness, but I found that this is actually superfluous,
since the intervals are all merged before it ever gets called. 
'''
def mergeOverlaps(intervalsCovered):
    retints = [intervalsCovered[0]] # start with the first interval
    for start, end in intervalsCovered:
        # if the current end is after the last intervals start and the current start
        # is before the end, we have either partial or full overlap
        if end >= retints[-1][0] and start <= retints[-1][1]:
            retints[-1][1] = max(end, retints[-1][1]) # expand the end of the interval
            retints[-1][0] = min(start, retints[-1][0]) # expand the beginning too
        else: # no overlap
            retints.append([start, end])
    return retints

if __name__=="__main__":
    # create the argument parser
    parser = argparse.ArgumentParser(description="Compute coverage statistics of a .bed file over refseq")
    parser.add_argument("reffile", help="The sorted reference file. Should be sorted according to chromosome, strand, then transcript start.")
    parser.add_argument("ref_format", help="The format of the reference file. Should be a one-line file that contains a header for the reference file (which should not contain any header) so that the script knows which column number refers to which column")
    parser.add_argument("bedfile", help="The .bed file that you want to compute overlap statistics for. Must be sorted according to chromosome, strand, then transcript start.")
    parser.add_argument("output", help="The name of the file you want to write the statistics to.")
    parser.add_argument("logfile", help="The name of the file you want to send the log to.")
    pargs = parser.parse_args()

    ## create the reference dict
    with open(pargs.ref_format, 'r') as formatfile:
        headerline = formatfile.readline().strip()
        headerdata = headerline.split("\t")
        ## from the UCSC browser, these column names will have prefixes of which reference
        ## table they came from, like mm9.knownCanonical. We don't care about these, so we
        ## look at the last entry split by periods
        ref_dict = {}
        for i in range(len(headerdata)):
            ref_dict[headerdata[i].split(".")[-1]] = i
                    
    compute_statistics(pargs.reffile, ref_dict, pargs.bedfile, pargs.output, pargs.logfile)

