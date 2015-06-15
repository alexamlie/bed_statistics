#!/usr/bin/python
## merge_reference_beds.py
## alex amlie-wolf 06-11-2015
## takes bed files and merges overlapping entries within each one
## for now, it is hard coded to use promoter, exon, intron, and intergenic files, but this can be
## easily extended (and the order doesn't actually matter)
## should be in the same directory as the chipseq.py file

import chipseq, sys

if __name__=="__main__":
    if len(sys.argv)==5:
        promoter_f = sys.argv[1]
        exon_f = sys.argv[2]
        intron_f = sys.argv[3]
        intergenic_f = sys.argv[4]
        
        chipseq.mergeBed(promoter_f)
        chipseq.mergeBed(exon_f)
        chipseq.mergeBed(intron_f)
        chipseq.mergeBed(intergenic_f)
    else:
        print "Usage: "+sys.argv[0]+" PROMOTER_FILE EXON_FILE INTRON_FILE INTERGENIC_FILE"
