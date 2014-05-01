'''
multicov_to_bedgraph.py
Alex Amlie-Wolf
A program to convert the output from a bedtools multicov file into a series of
bedGraph files, for visualization on the genome browser
'''

import sys

def convert_multicov(cov_file, outdir):
    cov = open(cov_file, 'r')
    header = cov.readline().strip().split("\t")
    ## the first 5 entries are chr, start, end, name, score, and strand
    exps = header[6:]
    ## make a dict of file descriptors
    exp_files = {}
    for i in exps:
        exp_files[i] = open(outdir+i+"_coverage.bedgraph", 'w')

    ## hard code different colors in the NLS experiment
    exp_colors = {"bg": "000,000,200", "nt":"200,000,000"}
        
    ## write the first line of the bedgraph format
    for exp in exp_files:
        if exp[:2] in exp_colors.keys():
            formatline = " ".join(["track type=bedGraph", "name="+exp, "visibility=full", "color="+exp_colors[exp[:2]]])
        else:
            formatline = " ".join(["track type=bedGraph", "name="+exp, "visibility=full", "color=000,000,200"])
        exp_files[exp].write(formatline+"\n")

    ## we want to keep track of which exons/introns we've seen already
    ## this will need to be optimized if i use this on a bigger file
    start_list = []
    end_list = []    
        
    ## now split the actual entries
    for entry in cov:
        data = entry.strip().split("\t")
        ## only write the data if we haven't seen it already
        if data[1] not in start_list and data[2] not in end_list:
            start_list.append(data[1])
            end_list.append(data[2])
            chrdata = "\t".join(data[:3])
            for i in range(len(exps)):
                ## write the chromosome data and then the counts (in the right order) 
                exp_files[exps[i]].write("\t".join([chrdata, data[i+6]])+"\n")
        
    [exp_files[i].close() for i in exp_files]
    cov.close()

if __name__=="__main__":
    if len(sys.argv) != 3:
        print "Usage: "+sys.argv[0]+" <multicov file> <output directory>"
    else:
        convert_multicov(sys.argv[1], sys.argv[2])
