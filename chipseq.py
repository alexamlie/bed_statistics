#!/usr/bin/python
# Contains useful functions for manipulating chIPseq data
# Greg Donahue, 06-16-2010
# ------------------------------------------------------------------------------
import sys
# ------------------------------------------------------------------------------
# GLOBALS
# Genome sizes
genome_sizes = { "hg18":3107677273, "hg19":3137161264,
                 "mm8":2664455088, "mm9":2725765481,
                 "sacCer1":12156302, "sacCer2":12156677 }

# ------------------------------------------------------------------------------
# CLASSES
# The BedReader class lets us easily spool through BED files
class BedReader:
    
    # Instance variables
    # handle is the file being read
    # delimiter is the record delimiter character
    # header is any header information collected from this file
    # spool holds the current line
    
    # Initialize new BedReaders
    def __init__(self, filename):
        self.handle = open(filename)
        self.delimiter = "\t" if self.usesTabs(filename) else " "
        self.header = dict()
        if self.hasHeader(filename): self.getHeader(self.handle.readline())
        self.spool = self.handle.readline()

    # Determine whether a BED file has a header
    def hasHeader(self, filename):
        f = open(filename); line = f.readline(); f.close()
        return line[0:5] == "track"

    # Determine whether a BED file uses tab or space delimiters
    def usesTabs(self, filename):
        f = open(filename); line = f.readline(); line = f.readline(); f.close()
        return line.count("\t") > 0

    # Get header information from the header string
    def getHeader(self, hstring):
        key, value, adding_to_key, in_quotes = "", "", True, False
        for i in range(6, len(hstring)-1):
            if hstring[i] == " ":
                if in_quotes: value += hstring[i]
                else:
                    self.header[key] = value
                    key, value, adding_to_key = "", "", True
            elif hstring[i] == "\"":
                if in_quotes:
                    self.header[key] = value+hstring[i]
                    key, value, adding_to_key, in_quotes = "", "", True, False
                else: value, in_quotes = value+hstring[i], True
            elif hstring[i] == "=": adding_to_key = False
            else:
                if adding_to_key: key += hstring[i]
                else: value += hstring[i]
        self.header[key] = value
        if "" in self.header.keys(): del self.header[""]

    # Return the next line
    def read(self):
        ret = self.format(self.spool)
        self.spool = self.handle.readline()
        return ret

    # Format a record string into a list of fields
    def format(self, line):
        ret, tokens = list(), line[0:-1].split(self.delimiter)
        try:
            ret.append(tokens[0])
            ret.append(int(tokens[1]))
            ret.append(int(tokens[2]))
            ret.append(tokens[3])
            ret.append(float(tokens[4]))
            ret.append(tokens[5])
            ret.append(int(tokens[6]))
            ret.append(int(tokens[7]))
            ret.append(tokens[8])
            ret.append(int(tokens[9]))
            ret.append([ int(bsize) for bsize in tokens[10].split(",") ])
            ret.append([ int(btart) for bstart in tokens[11].split(",") ])
        except Exception, e:
            try: ret.append(tokens[4])
            except Exception, e: ret.append(".")
        return ret

    # Does the file have any more data?
    def hasMoreData(self): return self.spool != ""
        
    # Close the file handle
    def close(self): self.handle.close()

    # Make a header string from the header dictionary
    def getHeaderString(self):
        ret = "track"
        for k in sorted(self.header.keys()): ret += " "+k+"="+self.header[k]
        return ret+"\n"

    # Make a record string from a record
    def getRecordString(self, record):
        ret = ""
        try:
            ret += record[0]
            ret += "\t"+str(record[1])
            ret += "\t"+str(record[2])
            ret += "\t"+record[3]
            ret += "\t"+str(record[4])
            ret += "\t"+record[5]
            ret += "\t"+str(record[6])
            ret += "\t"+str(record[7])
            ret += "\t"+record[8]
            ret += "\t"+str(record[9])
            ret += "\t"+str(record[10][0])
            for bsize in record[10][1:]: ret += ","+str(bsize)
            ret += "\t"+str(record[11][0])
            for bstart in record[11][1:]: ret += ","+str(bstart)
        except Exception, e: pass
        return ret+"\n"

# ------------------------------------------------------------------------------
# FUNCTIONS
# Convert a FASTA file downloaded from UCSC GB to one PWMSCAN can parse
def convertFasta(filename):
    fasta = getFasta(filename)
    f = open(filename[0:-3]+".converted.fa", 'w')
    for k in sorted(fasta.keys()):
        f.write(">"+k.split(" ")[1]+"\n"+fasta[k]+"\n")
    f.close()

# Get a FASTA object (dictionary) from a FASTA file
def getFasta(filename):
    ret = dict()
    f = open(filename)
    line = f.readline()
    while line != "":
        if line[0] == ">":
            header = line[1:-1]
            ret[header] = ""
        else: ret[header] += line[0:-1]
        line = f.readline()
    f.close()
    return ret

# Turn a GFF file generated by PWMSCAN into a BED file
def gffToBed(filename):
    gff = dict()
    f = open(filename); line = f.readline()
    while line != "":
        t = line[0:-1].split("\t")
        chromosome = t[0].split("=")[1].split(":")[0]
        start = int(t[0].split("=")[1].split(":")[1].split("-")[0])
        stop = start+int(t[4])
        start += int(t[3])
        try: gff[chromosome].append((start,stop))
        except Exception, e: gff[chromosome] = [ (start,stop) ]
        line = f.readline()
    f.close()
    f = open(filename[0:-3]+"bed", 'w')
    for chromosome in gff.keys():
        for record in gff[chromosome]:
            f.write(chromosome+"\t"+str(record[0])+"\t"+str(record[1])+"\n")
    f.close()

# Rewrite bowtie output as a BED file
def bowtieToBed(filename, name="ChIPseq", description="ChIPseq Raw Reads",
                color="10,10,10"):
    f = open(filename); line = f.readline()
    g = open(filename[0:-3]+"bed", 'w')
    #g.write("track name=\""+name+"\" description=\""+
    #        description+"\" color="+color+" visibility=full\n")
    while line != "":
        t = line[0:-1].split("\t")
        g.write(t[2]+"\t"+t[3]+"\t"+str(int(t[3])+len(t[4]))+"\t"+t[0]+
                "\t0\t"+t[1]+"\n")
        line = f.readline()
    f.close()
    g.close()

# Get global statistics on a chIPseq file
def getGlobalStatistics(filename, genome_size):
    br = BedReader(filename)
    count, average = 0.0, 0.0
    while br.hasMoreData():
        record = br.read()
        count += 1
        average += record[2]-record[1]
    br.close()
    average /= count
    return { "Count":count,
             "Average Tag Size":average,
             "RPKM Coefficient":(1000*1000000)/(average*count),
             "Coverage":count*average/genome_size }

# Clean a BED file for preprocessing by BEDtools, etc
def cleanBed(filename):
    br = BedReader(filename); g = open(filename[0:-3]+"clean.bed", 'w')
    while br.hasMoreData(): g.write(br.getRecordString(br.read()))
    br.close(); g.close()

# Get overlaps between two BED files
def getOverlaps(bedfile_a, bedfile_b):
    ret = dict()
    primary_loci = loadLocusDictionary(bedfile_a)
    secondary_loci = loadLocusDictionary(bedfile_b)
    for chromosome in primary_loci.keys():
        if not chromosome in secondary_loci.keys(): continue
        for primary in primary_loci[chromosome]:
            for secondary in secondary_loci[chromosome]:
                if areOverlapping(primary, secondary): ret[tuple(primary)] = 0
    return sorted(ret.keys())

# Determine whether two loci (BED records) are overlapping
def areOverlapping(a, b, symmetric_overlap=False):
    if a[1] >= b[1] and a[1] <= b[2]: return True
    if a[2] >= b[1] and a[2] <= b[2]: return True
    if b[1] >= a[1] and b[1] <= a[2]: return True
    if b[2] >= a[1] and b[2] <= a[2]: return True
    return False

# Load a locus dictionary of the form CHROMOSOME->[ RECORD_1, ..., RECORD_N ]
def loadLocusDictionary(filename):
    ret = dict()
    br = BedReader(filename)
    while br.hasMoreData():
        record = br.read()
        try: ret[record[0]].append(record)
        except Exception, e: ret[record[0]] = [ record ]
    br.close()
    return ret

# Count the number of entries in a locus dictionary
def countLocusDictionary(ldict):
    ret = 0
    for chromosome in ldict.keys(): ret += len(ldict[chromosome])
    return ret

# Load an empty bin dictionary from a BED file
def loadEmptyBinDictionary(filename, bin_size):
    br = BedReader(filename)
    ret = dict()
    while br.hasMoreData():
        record = br.read()
        i = int(record[1]/bin_size)*bin_size
        while i < int(record[2]/bin_size)*bin_size:
            try: ret[record[0]][i] = 0
            except Exception, e: ret[record[0]] = { i:0 }
            i += bin_size
    br.close()
    return ret

# Add reads from a BED file to a bin dictionary
def loadBinDictionary(filename, bins, bin_size):
    br = BedReader(filename)
    while br.hasMoreData():
        record = br.read()
        if len(record) > 5: start = record[1] if record[5] == "+" else record[2]
        else: start = record[1]
        try: bins[record[0]][int(start/bin_size)*bin_size] += 1
        except Exception, e: pass
    br.close()

# Normalize a bin dictionary by its RPKM
def normalizeBinDictionary(bins, rpkm):
    for chromosome in bins.keys():
        for bin in bins[chromosome].keys():
            bins[chromosome][bin] *= rpkm

# Subtract one bin dictionary from another
def subtractBinDictionaries(bins, background):
    for chromosome in background.keys():
        if not chromosome in bins.keys(): continue
        for bin in background[chromosome].keys():
            try: bins[chromosome][bin] -= background[chromosome][bin]
            except Exception, e:
                bins[chromosome][bin] = -1*background[chromosome][bin]

# Merge overlapping regions in a BED file and write another BED
def mergeBed(filename):
    regions = dict()
    br = BedReader(filename)
    while br.hasMoreData():
        r = br.read()
        try: regions[r[0]].append((r[0],r[1],r[2]))
        except Exception, e: regions[r[0]] = [ (r[0],r[1],r[2]) ]
    br.close()
    merged = dict()
    for chromosome in regions.keys():
        merged[chromosome] = list()
        loci = sorted(regions[chromosome])
        previous = loci[0]
        for locus in loci[1:]:
            if areOverlapping(previous, locus):
                previous = (previous[0],previous[1],locus[2])
            else:
                merged[chromosome].append(previous)
                previous = locus
        merged[chromosome].append(locus)
    f = open(filename[:-3]+"merged.bed", 'w')
    for chromosome in sorted(merged.keys()):
        for locus in merged[chromosome]:
            f.write(chromosome+"\t"+str(locus[1])+"\t"+str(locus[2])+"\n")
    f.close()

# Prune a BED file to remove all entries with out-of-bounds errors
def pruneBED(filename, size_file):
    sizes = loadSizes(size_file)
    br = BedReader(filename)
    f = open(filename[:-3]+"pruned.bed", 'w')
    while br.hasMoreData():
        r = br.read()
        if r[1] >= 0 and r[2] <= sizes[r[0]]: f.write(br.getRecordString(r))
    br.close()
    f.close()

# Prune a BedGraph file to remove all entries with out-of-bounds errors
def pruneBGR(filename, size_file):
    sizes = loadSizes(size_file)
    f = open(filename); line = f.readline()
    g = open(filename[:-3]+"pruned.bgr", 'w')
    while line != "":
        if len(line) > 0 and not line[0] in [ "b", "t", "#" ]:
            t = line[:-1].split("\t")
            if int(t[1]) >= 0 and int(t[2]) <= sizes[t[0]]: g.write(line)
        else: g.write(line)
        line = f.readline()
    f.close()
    g.close()

# Load a size file from UCSC Genome Browser
def loadSizes(filename):
    ret = dict()
    f = open(filename); lines = f.readlines(); f.close()
    for line in lines:
        t = line[:-1].split("\t")
        ret[t[0]] = int(t[1])
    return ret

# Load a dictionary of the form (chromosome,start,stop)->score from a bedGraph
def loadBedGraph(filename):
    f = open(filename); lines = f.readlines(); f.close()
    i = 0
    while i < len(lines):
        if lines[i][0] in [ "t", "b", "#" ]: i += 1
        else: break
    ret = dict()
    for line in lines[i:]:
        t = line[:-1].split("\t")
        ret[(t[0],int(t[1]),int(t[2]))] = float(t[3])
    return ret

# ------------------------------------------------------------------------------
# The following code is executed upon command-line invocation
if __name__ == "__main__": main(sys.argv)

# ------------------------------------------------------------------------------
# EOF
