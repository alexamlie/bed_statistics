#!/usr/bin/python

## merge_bed_files.py
## alex amlie-wolf 07-15-2015
## takes bed files and merges overlapping entries within each one
## unlike merge_reference_beds.py, this takes an arbitrary number of bed files and just
## merges all of them (does not check whether they are valid or anything)
## should be in the same directory as the chipseq.py file

import chipseq, sys

if __name__=="__main__":
    if len(sys.argv) > 2:
        for i in range(1, len(sys.argv)):
            chipseq.mergeBed(sys.argv[i])
    else:
        print "Usage: "+sys.argv[0]+" <bed file 1> <bed file 2> ...."
        print "(Takes an arbitrary number of files)"        
