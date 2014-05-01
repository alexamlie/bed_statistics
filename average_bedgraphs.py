'''
average_bedgraphs.py
Alex Amlie-Wolf, Spring 2014
A program to take in sets of bedgraph files and compute the average coverage of conditions
Note that the output file must be given first for the nargs="*" to work
Assumes that the bedgraph files are all identical in each condition
'''

import argparse

def compute_averages(cond1_files, cond2_files, c1_output, c2_output, c1name, c2name):
    c1out = open(c1_output, 'w')
    c2out = open(c2_output, 'w')
    cond1 = {}
    cond2 = {}
    for i in range(len(cond1_files)):
        cond1[cond1_files[i]] = open(cond1_files[i], 'r')
    for i in range(len(cond2_files)):
        cond2[cond2_files[i]] = open(cond2_files[i], 'r')

    ## do condition 1 first
    ## first read through the headers
    for f in cond1_files:
        headerdata = cond1[f].readline().strip().split(" ")
        color = headerdata[-1].split("=")[1]        
                
    c1out.write(" ".join(["track type=bedGraph", "name="+c1name+"_average", "visibility=full", "color="+color])+"\n") 
    
    c1_done = False
    while not c1_done:
        coverage_sum = 0.0
        for f in cond1_files:
            line = cond1[f].readline().strip()
            if line:
                data = line.split("\t")
                coverage_sum += int(data[3])
                chrdata = "\t".join(data[:3])
            else:
                c1_done = True
                break

        if not c1_done:            
            avg_coverage = coverage_sum / float(len(cond1_files))        
            c1out.write("\t".join([chrdata, str(avg_coverage)])+"\n")
    
    for f in cond2_files:
        headerdata = cond2[f].readline().strip().split(" ")
        color = headerdata[-1].split("=")[1]
                
    c2out.write(" ".join(["track type=bedGraph", "name="+c2name+"_average", "visibility=full", "color="+color])+"\n") 
    
    c2_done = False
    while not c2_done:
        coverage_sum = 0.0
        for f in cond2_files:
            line = cond2[f].readline().strip()
            if line:
                data = line.split("\t")
                coverage_sum += int(data[3])
                chrdata = "\t".join(data[:3])
            else:
                c2_done = True

        if not c2_done:
            avg_coverage = coverage_sum / float(len(cond2_files))
            c2out.write("\t".join([chrdata, str(avg_coverage)])+"\n")
                
    c1out.close()
    c2out.close()
    [cond1[i].close() for i in cond1]
    [cond2[i].close() for i in cond2]

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Get averages of bedgraph files")
    parser.add_argument("--c1", nargs="*", help="The bedgraph files in condition 1")
    parser.add_argument("--c2", nargs="*", help="The bedgraph files in condition 2")
    parser.add_argument("--c1name", help="The name of condition 1")
    parser.add_argument("--c2name", help="The name of condition 2")
    parser.add_argument("c1_output", help="The output file for condition 1")
    parser.add_argument("c2_output", help="The output file for condition 2")   
    pargs = parser.parse_args()
    compute_averages(pargs.c1, pargs.c2, pargs.c1_output, pargs.c2_output, pargs.c1name, pargs.c2name)
