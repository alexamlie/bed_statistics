"""
entrywise_bed_coverage.py
Alex Amlie-Wolf, 6/10/15

A program to compute coverage statistics on each entry of a .bed file given bed files for
promoters, exons, introns, and repeat regions, as well as optional files for either full UTR
regions or exons and introns. All .bed files MUST be sorted by chromosome and strand, otherwise
we will miss peaks. It outputs one file containing each entry of the input bed file with its
overlap with each type of element, and another file containing summary statistics for how much
overlap there was for each element

For now, takes in bed files with only three columns (because bedtools complement removes strand
information). To do stranded analysis, you have to run the script
separately using element information for each strand.

Requires Python >= 2.7
"""

import argparse, sys, os, time, gzip, re

def compute_coverage(promoter_f, exon_f, intron_f, repeat_f, input_f, entrywise_out_f, summary_out_f):
    start = time.clock()
    """
    the main function to compute the coverage of the entries in the bed file.
    """
    # define a convenience dictionary to make subsetting clearer
    bed_coords = {} 
    bed_coords["chrom"] = 0
    bed_coords["start"] = 1
    bed_coords["end"] = 2    
    bed_coords["strand"] = 5

    if input_f[-2:]=='gz':
        input_beds = gzip.open(input_f, 'rb')
    else:
        input_beds = open(input_f, 'r')
    
    with open(promoter_f, 'r') as promoters, open(exon_f, 'r') as exons, open(intron_f, 'r') as introns, open(repeat_f, 'r') as repeats, open(entrywise_out_f, 'w') as entry_out, open(summary_out_f, 'w') as summary_out:
        # loop through input bed file, and store lists of potential peaks from other file types
        # (we don't want to stop considering a peak until the start site that we are looking at is
        # past its end site)

        ## write the output file headers
        # for the individual entries, write the original entry along with its amount and percent of
        # overlap with each element
        ## might want to add in 'other' here, although it's implied 
        entry_out.write("\t".join(['chr', 'start', 'end', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')
        # for the summary, we write each chromosomes entry as well as genomewide
        summary_out.write('\t'.join(['partition', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')

        ## we store lists of each type of element
        cur_promoters = []
        cur_exons = []
        cur_introns = []
        cur_repeats = []

        ## start reading in these genomic elements
        cur_promoter = promoters.readline().strip().split('\t')
        cur_exon = exons.readline().strip().split('\t')
        cur_intron = introns.readline().strip().split('\t')
        cur_repeat = repeats.readline().strip().split('\t')

        ## also initialize our tracking variables for the summary statistics:
        total_entry_bp = 0.0 # the total number of base pairs covered by the input bed
        total_promoter_bp = 0.0 # the total number of promoter BPs overlapped
        total_exon_bp = 0.0 # total number of exon BPs overlapped
        total_intron_bp = 0.0 # total number of intron BPs overlapped
        total_repeat_bp = 0.0 # total number of repeat BPs overlapped

        # for each chromosome, also store similar info
        this_chr_entry_bp = 0.0
        this_chr_promoter_bp = 0.0 
        this_chr_exon_bp = 0.0 
        this_chr_intron_bp = 0.0 
        this_chr_repeat_bp = 0.0 

        # to track which chromosome we're on
        this_chr = ''

        ## need to make sure that if the entry's chromosome is not in our reference file,
        ## we just wait..do this by checking if this_chr is in the form of
        ## chr[0-9]_* (i.e. it's got some random crap after it)
        ## i am assuming that the reference files only contain canonical chromosomes
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        
                
        entry_ctr = 0        
        for entry in input_beds:
            entry_ctr += 1
            if entry_ctr % 1000==0:
                print entry_ctr
            entry_data = entry.strip().split('\t')
            entry_start = int(entry_data[bed_coords['start']])
            entry_end = int(entry_data[bed_coords['end']])
            this_length = int(entry_data[bed_coords['end']]) - entry_start
            total_entry_bp += this_length
            if entry_data[bed_coords['chrom']] == this_chr:
                this_chr_entry_bp += this_length
            else:
                if this_chr: # if we've read anything in yet:
                    summary_out.write('\t'.join([this_chr, str(this_chr_promoter_bp), str(this_chr_promoter_bp/this_chr_entry_bp), str(this_chr_exon_bp), str(this_chr_exon_bp/this_chr_entry_bp), str(this_chr_intron_bp), str(this_chr_intron_bp/this_chr_entry_bp), str(this_chr_repeat_bp), str(this_chr_repeat_bp/this_chr_entry_bp)])+'\n')
                ## reset the genomic element lists and read through their files until they're
                ## on the right chromosome again
                this_chr = entry_data[bed_coords['chrom']]
                cur_promoters = []
                cur_exons = []
                cur_introns = []
                cur_repeats = []
                ## reset the chromosome tracking variables
                this_chr_entry_bp = float(this_length)
                this_chr_promoter_bp = 0.0 
                this_chr_exon_bp = 0.0 
                this_chr_intron_bp = 0.0 
                this_chr_repeat_bp = 0.0                                 

                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    continue
                
                while len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']] != this_chr:
                    cur_promoter = promoters.readline().strip().split('\t')
                while len(cur_exon) > 1 and cur_exon[bed_coords['chrom']] != this_chr:
                    cur_exon = exons.readline().strip().split('\t')
                while len(cur_intron) > 1 and cur_intron[bed_coords['chrom']] != this_chr:
                    cur_intron = introns.readline().strip().split('\t')
                while len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']] != this_chr:
                    cur_repeat = repeats.readline().strip().split('\t')
                print 'Parsing chromosome '+this_chr
            ## add all possibly overlapping entries for each type of element
            ## must have same chromosome and start before the end of the entry
            promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
            while promoterbool:
                cur_promoters.append(cur_promoter)
                cur_promoter = promoters.readline().strip().split('\t')
                promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
                
            exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end
            while exonbool:
                cur_exons.append(cur_exon)
                cur_exon = exons.readline().strip().split('\t')
                exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end 

            intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 
            while intronbool:
                cur_introns.append(cur_intron)
                cur_intron = introns.readline().strip().split('\t')
                intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 

            repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 
            while repeatbool:
                cur_repeats.append(cur_repeat)
                cur_repeat = repeats.readline().strip().split('\t')
                repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 

            ## remove entries that are before the current entry here:
            for i in range(len(cur_promoters)-1, -1, -1):
                if cur_promoters[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_promoters[i][bed_coords['end']]) < entry_start:
                    cur_promoters.pop(i)

            for i in range(len(cur_exons)-1, -1, -1):
                if cur_exons[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_exons[i][bed_coords['end']]) < entry_start:
                    cur_exons.pop(i)

            for i in range(len(cur_introns)-1, -1, -1):
                if cur_introns[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_introns[i][bed_coords['end']]) < entry_start:
                    cur_introns.pop(i)

            for i in range(len(cur_repeats)-1, -1, -1):
                if cur_repeats[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_repeats[i][bed_coords['end']]) < entry_start:
                    cur_repeats.pop(i)                                        
                                    
            ## now compare the current entry to all these peaks
            ## here i am relying on the correct merging of the bed entries from the element files
            ## because i dont ever re-merge anything
            this_promoter_bp = 0.0
            this_exon_bp = 0.0
            this_intron_bp = 0.0
            this_repeat_bp = 0.0
            for cp in cur_promoters:
                # if int(cp[bed_coords['start']]) > int(entry_data[bed_coords['end']]):
                #     break
                this_promoter_bp += max(0, min(int(cp[bed_coords['end']]), entry_end) - max(int(cp[bed_coords['start']]), entry_start))
            for ce in cur_exons:
                # if int(ce[bed_coords['start']]) > int(entry_data[bed_coords['end']]):
                #     break
                this_exon_bp += max(0, min(int(ce[bed_coords['end']]), entry_end) - max(int(ce[bed_coords['start']]), entry_start))
            for ci in cur_introns:
                # if int(ci[bed_coords['start']]) > int(entry_data[bed_coords['end']]):
                #     break                
                this_intron_bp += max(0, min(int(ci[bed_coords['end']]), entry_end) - max(int(ci[bed_coords['start']]), entry_start))
            for cr in cur_repeats:
                # if int(cr[bed_coords['start']]) > int(entry_data[bed_coords['end']]):
                #     break                
                this_repeat_bp += max(0, min(int(cr[bed_coords['end']]), entry_end) - max(int(cr[bed_coords['start']]), entry_start))
            ## add these overlaps to the summary trackers
            this_chr_promoter_bp += this_promoter_bp
            this_chr_exon_bp += this_exon_bp
            this_chr_intron_bp += this_intron_bp
            this_chr_repeat_bp += this_repeat_bp
            total_promoter_bp += this_promoter_bp
            total_exon_bp += this_exon_bp
            total_intron_bp += this_intron_bp
            total_repeat_bp += this_repeat_bp
            
            entry_out.write('\t'.join([entry_data[bed_coords['chrom']], entry_data[bed_coords['start']], entry_data[bed_coords['end']], str(this_promoter_bp), str(this_promoter_bp / this_length), str(this_exon_bp), str(this_exon_bp / this_length), str(this_intron_bp), str(this_intron_bp / this_length), str(this_repeat_bp), str(this_repeat_bp / this_length)])+'\n')

        summary_out.write('\t'.join(['genomewide', str(total_promoter_bp), str(total_promoter_bp/total_entry_bp), str(total_exon_bp), str(total_exon_bp/total_entry_bp), str(total_intron_bp), str(total_intron_bp/total_entry_bp), str(total_repeat_bp), str(total_repeat_bp/total_entry_bp)])+'\n')

    input_beds.close()
        
    end = time.clock()
    length = end - start
    print "Analysis complete, time: ", length

## another function to calculate coverage over the UTRs (given one file for each UTR):
def compute_full_utr_coverage(fp_utr_f, tp_utr_f, promoter_f, exon_f, intron_f, repeat_f, input_f, entrywise_out_f, summary_out_f):
    start = time.clock()
    """
    the main function to compute the coverage of the entries in the bed file.
    """
    # define a convenience dictionary to make subsetting clearer
    bed_coords = {} 
    bed_coords["chrom"] = 0
    bed_coords["start"] = 1
    bed_coords["end"] = 2    
    bed_coords["strand"] = 5

    if input_f[-2:]=='gz':
        input_beds = gzip.open(input_f, 'rb')
    else:
        input_beds = open(input_f, 'r')
    
    with open(fp_utr_f, 'r') as fp_utrs, open(tp_utr_f, 'r') as tp_utrs, open(promoter_f, 'r') as promoters, open(exon_f, 'r') as exons, open(intron_f, 'r') as introns, open(repeat_f, 'r') as repeats, open(entrywise_out_f, 'w') as entry_out, open(summary_out_f, 'w') as summary_out:
        # loop through input bed file, and store lists of potential peaks from other file types
        # (we don't want to stop considering a peak until the start site that we are looking at is
        # past its end site)

        ## write the output file headers
        # for the individual entries, write the original entry along with its amount and percent of
        # overlap with each element
        entry_out.write("\t".join(['chr', 'start', 'end', 'fp_utr_bp', 'fp_utr_pct', 'tp_utr_bp', 'tp_utr_pct', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')
        # for the summary, we write each chromosomes entry as well as genomewide
        summary_out.write('\t'.join(['partition', 'fp_utr_bp', 'fp_utr_pct', 'tp_utr_bp', 'tp_utr_pct', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')

        ## we store lists of each type of element
        cur_fp_utrs = []
        cur_tp_utrs = []
        cur_promoters = []
        cur_exons = []
        cur_introns = []
        cur_repeats = []

        ## start reading in these genomic elements
        cur_fp_utr = fp_utrs.readline().strip().split('\t')
        cur_tp_utr = tp_utrs.readline().strip().split('\t')                
        cur_promoter = promoters.readline().strip().split('\t')
        cur_exon = exons.readline().strip().split('\t')
        cur_intron = introns.readline().strip().split('\t')
        cur_repeat = repeats.readline().strip().split('\t')

        ## also initialize our tracking variables for the summary statistics:
        total_entry_bp = 0.0 # the total number of base pairs covered by the input bed
        total_fp_utr_bp = 0.0 # total number of 5' UTR bps
        total_tp_utr_bp = 0.0 # total number of 3' UTR bps        
        total_promoter_bp = 0.0 # the total number of promoter BPs overlapped
        total_exon_bp = 0.0 # total number of exon BPs overlapped
        total_intron_bp = 0.0 # total number of intron BPs overlapped
        total_repeat_bp = 0.0 # total number of repeat BPs overlapped

        # for each chromosome, also store similar info
        this_chr_entry_bp = 0.0
        this_chr_fp_utr_bp = 0.0
        this_chr_tp_utr_bp = 0.0
        this_chr_promoter_bp = 0.0 
        this_chr_exon_bp = 0.0 
        this_chr_intron_bp = 0.0 
        this_chr_repeat_bp = 0.0 

        # to track which chromosome we're on
        this_chr = ''

        ## need to make sure that if the entry's chromosome is not in our reference file,
        ## we just wait..do this by checking if this_chr is in the form of
        ## chr[0-9]_* (i.e. it's got some random crap after it)
        ## i am assuming that the reference files only contain canonical chromosomes
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        
        
        entry_ctr = 0        
        for entry in input_beds:
            entry_ctr += 1
            if entry_ctr % 1000==0:
                print entry_ctr
            entry_data = entry.strip().split('\t')
            entry_start = int(entry_data[bed_coords['start']])
            entry_end = int(entry_data[bed_coords['end']])
            this_length = int(entry_data[bed_coords['end']]) - entry_start
            total_entry_bp += this_length
            if entry_data[bed_coords['chrom']] == this_chr:
                this_chr_entry_bp += this_length
            else:
                if this_chr: # if we've read anything in yet:
                    summary_out.write('\t'.join([this_chr, str(this_chr_fp_utr_bp), str(this_chr_fp_utr_bp/this_chr_entry_bp), str(this_chr_tp_utr_bp), str(this_chr_tp_utr_bp/this_chr_entry_bp), str(this_chr_promoter_bp), str(this_chr_promoter_bp/this_chr_entry_bp), str(this_chr_exon_bp), str(this_chr_exon_bp/this_chr_entry_bp), str(this_chr_intron_bp), str(this_chr_intron_bp/this_chr_entry_bp), str(this_chr_repeat_bp), str(this_chr_repeat_bp/this_chr_entry_bp)])+'\n')
                ## reset the genomic element lists and read through their files until they're
                ## on the right chromosome again
                this_chr = entry_data[bed_coords['chrom']]
                cur_fp_utrs = []
                cur_tp_utrs = []                
                cur_promoters = []
                cur_exons = []
                cur_introns = []
                cur_repeats = []

                ## reset the chromosome tracking variables
                this_chr_entry_bp = float(this_length)
                this_chr_fp_utr_bp = 0.0
                this_chr_tp_utr_bp = 0.0                
                this_chr_promoter_bp = 0.0 
                this_chr_exon_bp = 0.0 
                this_chr_intron_bp = 0.0 
                this_chr_repeat_bp = 0.0                                 

                ## if we have a match to one of the weird chromosomes, just skip to the next entry
                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    continue
                
                while len(cur_fp_utr) > 1 and cur_fp_utr[bed_coords['chrom']] != this_chr:
                    cur_fp_utr = fp_utrs.readline().strip().split('\t')
                while len(cur_tp_utr) > 1 and cur_tp_utr[bed_coords['chrom']] != this_chr:
                    cur_tp_utr = tp_utrs.readline().strip().split('\t')                    
                while len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']] != this_chr:
                    cur_promoter = promoters.readline().strip().split('\t')
                while len(cur_exon) > 1 and cur_exon[bed_coords['chrom']] != this_chr:
                    cur_exon = exons.readline().strip().split('\t')
                while len(cur_intron) > 1 and cur_intron[bed_coords['chrom']] != this_chr:
                    cur_intron = introns.readline().strip().split('\t')
                while len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']] != this_chr:
                    cur_repeat = repeats.readline().strip().split('\t')
                print 'Parsing chromosome '+this_chr
            ## add all possibly overlapping entries for each type of element
            ## must have same chromosome and start before the end of the entry
            fp_utrbool = len(cur_fp_utr) > 1 and cur_fp_utr[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_utr[bed_coords['start']]) <= entry_end
            while fp_utrbool:
                cur_fp_utrs.append(cur_fp_utr)
                cur_fp_utr = fp_utrs.readline().strip().split('\t')
                fp_utrbool = len(cur_fp_utr) > 1 and cur_fp_utr[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_utr[bed_coords['start']]) <= entry_end

            tp_utrbool = len(cur_tp_utr) > 1 and cur_tp_utr[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_utr[bed_coords['start']]) <= entry_end
            while tp_utrbool:
                cur_tp_utrs.append(cur_tp_utr)
                cur_tp_utr = tp_utrs.readline().strip().split('\t')
                tp_utrbool = len(cur_tp_utr) > 1 and cur_tp_utr[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_utr[bed_coords['start']]) <= entry_end
                        
            promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
            while promoterbool:
                cur_promoters.append(cur_promoter)
                cur_promoter = promoters.readline().strip().split('\t')
                promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
                
            exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end
            while exonbool:
                cur_exons.append(cur_exon)
                cur_exon = exons.readline().strip().split('\t')
                exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end 

            intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 
            while intronbool:
                cur_introns.append(cur_intron)
                cur_intron = introns.readline().strip().split('\t')
                intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 

            repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 
            while repeatbool:
                cur_repeats.append(cur_repeat)
                cur_repeat = repeats.readline().strip().split('\t')
                repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 

            ## remove entries that are before the current entry here:
            for i in range(len(cur_fp_utrs)-1, -1, -1):
                if cur_fp_utrs[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_fp_utrs[i][bed_coords['end']]) < entry_start:
                    cur_fp_utrs.pop(i)

            for i in range(len(cur_tp_utrs)-1, -1, -1):
                if cur_tp_utrs[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_tp_utrs[i][bed_coords['end']]) < entry_start:
                    cur_tp_utrs.pop(i)
            
            for i in range(len(cur_promoters)-1, -1, -1):
                if cur_promoters[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_promoters[i][bed_coords['end']]) < entry_start:
                    cur_promoters.pop(i)

            for i in range(len(cur_exons)-1, -1, -1):
                if cur_exons[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_exons[i][bed_coords['end']]) < entry_start:
                    cur_exons.pop(i)

            for i in range(len(cur_introns)-1, -1, -1):
                if cur_introns[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_introns[i][bed_coords['end']]) < entry_start:
                    cur_introns.pop(i)

            for i in range(len(cur_repeats)-1, -1, -1):
                if cur_repeats[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_repeats[i][bed_coords['end']]) < entry_start:
                    cur_repeats.pop(i)                                        
                                    
            ## now compare the current entry to all these peaks
            ## here i am relying on the correct merging of the bed entries from the element files
            ## because i dont ever re-merge anything
            this_fp_utr_bp = 0.0
            this_tp_utr_bp = 0.0            
            this_promoter_bp = 0.0
            this_exon_bp = 0.0
            this_intron_bp = 0.0
            this_repeat_bp = 0.0
            for cfp in cur_fp_utrs:
                this_fp_utr_bp += max(0, min(int(cfp[bed_coords['end']]), entry_end) - max(int(cfp[bed_coords['start']]), entry_start))
            for ctp in cur_tp_utrs:
                this_tp_utr_bp += max(0, min(int(ctp[bed_coords['end']]), entry_end) - max(int(ctp[bed_coords['start']]), entry_start))                                
            for cp in cur_promoters:
                this_promoter_bp += max(0, min(int(cp[bed_coords['end']]), entry_end) - max(int(cp[bed_coords['start']]), entry_start))
            for ce in cur_exons:
                this_exon_bp += max(0, min(int(ce[bed_coords['end']]), entry_end) - max(int(ce[bed_coords['start']]), entry_start))
            for ci in cur_introns:
                this_intron_bp += max(0, min(int(ci[bed_coords['end']]), entry_end) - max(int(ci[bed_coords['start']]), entry_start))
            for cr in cur_repeats:
                this_repeat_bp += max(0, min(int(cr[bed_coords['end']]), entry_end) - max(int(cr[bed_coords['start']]), entry_start))
            ## add these overlaps to the summary trackers
            this_chr_fp_utr_bp += this_fp_utr_bp
            this_chr_tp_utr_bp += this_tp_utr_bp                        
            this_chr_promoter_bp += this_promoter_bp
            this_chr_exon_bp += this_exon_bp
            this_chr_intron_bp += this_intron_bp
            this_chr_repeat_bp += this_repeat_bp
            
            total_fp_utr_bp += this_fp_utr_bp
            total_tp_utr_bp += this_tp_utr_bp            
            total_promoter_bp += this_promoter_bp
            total_exon_bp += this_exon_bp
            total_intron_bp += this_intron_bp
            total_repeat_bp += this_repeat_bp
            
            entry_out.write('\t'.join([entry_data[bed_coords['chrom']], entry_data[bed_coords['start']], entry_data[bed_coords['end']], str(this_fp_utr_bp), str(this_fp_utr_bp / this_length), str(this_tp_utr_bp), str(this_tp_utr_bp / this_length), str(this_promoter_bp), str(this_promoter_bp / this_length), str(this_exon_bp), str(this_exon_bp / this_length), str(this_intron_bp), str(this_intron_bp / this_length), str(this_repeat_bp), str(this_repeat_bp / this_length)])+'\n')

        summary_out.write('\t'.join(['genomewide', str(total_fp_utr_bp), str(total_fp_utr_bp / total_entry_bp), str(total_tp_utr_bp), str(total_tp_utr_bp / total_entry_bp), str(total_promoter_bp), str(total_promoter_bp/total_entry_bp), str(total_exon_bp), str(total_exon_bp/total_entry_bp), str(total_intron_bp), str(total_intron_bp/total_entry_bp), str(total_repeat_bp), str(total_repeat_bp/total_entry_bp)])+'\n')

    input_beds.close()
        
    end = time.clock()
    length = end - start
    print "Analysis complete, time: ", length

## another function to calculate coverage over the UTRs
## takes two files for each UTR, exons and introns:
def compute_split_utr_coverage(fp_utr_exons_f, fp_utr_introns_f, tp_utr_exons_f, tp_utr_introns_f, promoter_f, exon_f, intron_f, repeat_f, input_f, entrywise_out_f, summary_out_f):
    start = time.clock()
    """
    the main function to compute the coverage of the entries in the bed file.
    """
    # define a convenience dictionary to make subsetting clearer
    bed_coords = {} 
    bed_coords["chrom"] = 0
    bed_coords["start"] = 1
    bed_coords["end"] = 2    
    bed_coords["strand"] = 5

    if input_f[-2:]=='gz':
        input_beds = gzip.open(input_f, 'rb')
    else:
        input_beds = open(input_f, 'r')
    
    with open(fp_utr_exons_f, 'r') as fp_exons, open(fp_utr_introns_f, 'r') as fp_introns, open(tp_utr_exons_f, 'r') as tp_exons, open(tp_utr_introns_f, 'r') as tp_introns, open(promoter_f, 'r') as promoters, open(exon_f, 'r') as exons, open(intron_f, 'r') as introns, open(repeat_f, 'r') as repeats, open(entrywise_out_f, 'w') as entry_out, open(summary_out_f, 'w') as summary_out:
        # loop through input bed file, and store lists of potential peaks from other file types
        # (we don't want to stop considering a peak until the start site that we are looking at is
        # past its end site)

        ## write the output file headers
        # for the individual entries, write the original entry along with its amount and percent of
        # overlap with each element
        entry_out.write("\t".join(['chr', 'start', 'end', 'fp_utr_exon_bp', 'fp_utr_exon_pct', 'fp_utr_intron_bp', 'fp_utr_intron_pct', 'tp_utr_exon_bp', 'tp_utr_exon_pct', 'tp_utr_intron_bp', 'tp_utr_intron_pct', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')
        # for the summary, we write each chromosomes entry as well as genomewide
        summary_out.write('\t'.join(['partition', 'fp_utr_exon_bp', 'fp_utr_exon_pct', 'fp_utr_intron_bp', 'fp_utr_intron_pct', 'tp_utr_exon_bp', 'tp_utr_exon_pct', 'tp_utr_intron_bp', 'tp_utr_intron_pct', 'promoter_bp', 'promoter_pct', 'exon_bp', 'exon_pct', 'intron_bp', 'intron_pct', 'repeat_bp', 'repeat_pct'])+'\n')

        ## we store lists of each type of element
        cur_fp_exons = []
        cur_fp_introns = []
        cur_tp_exons = []
        cur_tp_introns = []        
        cur_promoters = []
        cur_exons = []
        cur_introns = []
        cur_repeats = []

        ## start reading in these genomic elements
        cur_fp_exon = fp_exons.readline().strip().split('\t')
        cur_fp_intron = fp_introns.readline().strip().split('\t')        
        cur_tp_exon = tp_exons.readline().strip().split('\t')
        cur_tp_intron = tp_introns.readline().strip().split('\t')        
        cur_promoter = promoters.readline().strip().split('\t')
        cur_exon = exons.readline().strip().split('\t')
        cur_intron = introns.readline().strip().split('\t')
        cur_repeat = repeats.readline().strip().split('\t')

        ## also initialize our tracking variables for the summary statistics:
        total_entry_bp = 0.0 # the total number of base pairs covered by the input bed
        total_fp_exon_bp = 0.0 # total number of 5' UTR exon bps
        total_fp_intron_bp = 0.0 # total number of 5' UTR intron bps        
        total_tp_exon_bp = 0.0 # total number of 3' UTR exon bps
        total_tp_intron_bp = 0.0 # total number of 3' UTR intron bps        
        total_promoter_bp = 0.0 # the total number of promoter BPs overlapped
        total_exon_bp = 0.0 # total number of exon BPs overlapped
        total_intron_bp = 0.0 # total number of intron BPs overlapped
        total_repeat_bp = 0.0 # total number of repeat BPs overlapped

        # for each chromosome, also store similar info
        this_chr_entry_bp = 0.0
        this_chr_fp_exon_bp = 0.0
        this_chr_fp_intron_bp = 0.0
        this_chr_tp_exon_bp = 0.0
        this_chr_tp_intron_bp = 0.0
        this_chr_promoter_bp = 0.0 
        this_chr_exon_bp = 0.0 
        this_chr_intron_bp = 0.0 
        this_chr_repeat_bp = 0.0 

        # to track which chromosome we're on
        this_chr = ''

        ## need to make sure that if the entry's chromosome is not in our reference file,
        ## we just wait..do this by checking if this_chr is in the form of
        ## chr[0-9]_* (i.e. it's got some random crap after it)
        ## i am assuming that the reference files only contain canonical chromosomes
        ## i also manually check against chrM later
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        
        
        entry_ctr = 0        
        for entry in input_beds:
            entry_ctr += 1
            if entry_ctr % 1000==0:
                print entry_ctr
            entry_data = entry.strip().split('\t')
            entry_start = int(entry_data[bed_coords['start']])
            entry_end = int(entry_data[bed_coords['end']])
            this_length = int(entry_data[bed_coords['end']]) - entry_start
            total_entry_bp += this_length
            if entry_data[bed_coords['chrom']] == this_chr:
                this_chr_entry_bp += this_length
            else:
                if this_chr: # if we've read anything in yet:
                    summary_out.write('\t'.join([this_chr, str(this_chr_fp_exon_bp), str(this_chr_fp_exon_bp/this_chr_entry_bp), str(this_chr_fp_intron_bp), str(this_chr_fp_intron_bp/this_chr_entry_bp), str(this_chr_tp_exon_bp), str(this_chr_tp_exon_bp/this_chr_entry_bp), str(this_chr_tp_intron_bp), str(this_chr_tp_intron_bp/this_chr_entry_bp), str(this_chr_promoter_bp), str(this_chr_promoter_bp/this_chr_entry_bp), str(this_chr_exon_bp), str(this_chr_exon_bp/this_chr_entry_bp), str(this_chr_intron_bp), str(this_chr_intron_bp/this_chr_entry_bp), str(this_chr_repeat_bp), str(this_chr_repeat_bp/this_chr_entry_bp)])+'\n')
                ## reset the genomic element lists and read through their files until they're
                ## on the right chromosome again
                this_chr = entry_data[bed_coords['chrom']]
                cur_fp_exons = []
                cur_fp_introns = []
                cur_tp_exons = []
                cur_tp_introns = []                                
                cur_promoters = []
                cur_exons = []
                cur_introns = []
                cur_repeats = []

                ## reset the chromosome tracking variables
                this_chr_entry_bp = float(this_length)
                this_chr_fp_exon_bp = 0.0
                this_chr_fp_intron_bp = 0.0
                this_chr_tp_exon_bp = 0.0
                this_chr_tp_intron_bp = 0.0                
                this_chr_promoter_bp = 0.0 
                this_chr_exon_bp = 0.0 
                this_chr_intron_bp = 0.0 
                this_chr_repeat_bp = 0.0                                 

                ## if we have a match to one of the weird chromosomes, just skip to the next entry
                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    continue
                
                while len(cur_fp_exon) > 1 and cur_fp_exon[bed_coords['chrom']] != this_chr:
                    cur_fp_exon = fp_exons.readline().strip().split('\t')
                while len(cur_fp_intron) > 1 and cur_fp_intron[bed_coords['chrom']] != this_chr:
                    cur_fp_intron = fp_introns.readline().strip().split('\t')
                while len(cur_tp_exon) > 1 and cur_tp_exon[bed_coords['chrom']] != this_chr:
                    cur_tp_exon = tp_exons.readline().strip().split('\t')
                while len(cur_tp_intron) > 1 and cur_tp_intron[bed_coords['chrom']] != this_chr:
                    cur_tp_intron = tp_introns.readline().strip().split('\t')
                while len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']] != this_chr:
                    cur_promoter = promoters.readline().strip().split('\t')
                while len(cur_exon) > 1 and cur_exon[bed_coords['chrom']] != this_chr:
                    cur_exon = exons.readline().strip().split('\t')
                while len(cur_intron) > 1 and cur_intron[bed_coords['chrom']] != this_chr:
                    cur_intron = introns.readline().strip().split('\t')
                while len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']] != this_chr:
                    cur_repeat = repeats.readline().strip().split('\t')
                print 'Parsing chromosome '+this_chr
                
            ## add all possibly overlapping entries for each type of element
            ## must have same chromosome and start before the end of the entry
            fp_exonbool = len(cur_fp_exon) > 1 and cur_fp_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_exon[bed_coords['start']]) <= entry_end
            while fp_exonbool:
                cur_fp_exons.append(cur_fp_exon)
                cur_fp_exon = fp_exons.readline().strip().split('\t')
                fp_exonbool = len(cur_fp_exon) > 1 and cur_fp_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_exon[bed_coords['start']]) <= entry_end
                
            fp_intronbool = len(cur_fp_intron) > 1 and cur_fp_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_intron[bed_coords['start']]) <= entry_end            
            while fp_intronbool:
                cur_fp_introns.append(cur_fp_intron)
                cur_fp_intron = fp_introns.readline().strip().split('\t')
                fp_intronbool = len(cur_fp_intron) > 1 and cur_fp_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_fp_intron[bed_coords['start']]) <= entry_end
                
            tp_exonbool = len(cur_tp_exon) > 1 and cur_tp_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_exon[bed_coords['start']]) <= entry_end            
            while tp_exonbool:
                cur_tp_exons.append(cur_tp_exon)
                cur_tp_exon = tp_exons.readline().strip().split('\t')
                tp_exonbool = len(cur_tp_exon) > 1 and cur_tp_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_exon[bed_coords['start']]) <= entry_end
                
            tp_intronbool = len(cur_tp_intron) > 1 and cur_tp_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_intron[bed_coords['start']]) <= entry_end            
            while tp_intronbool:
                cur_tp_introns.append(cur_tp_intron)
                cur_tp_intron = tp_introns.readline().strip().split('\t')
                tp_intronbool = len(cur_tp_intron) > 1 and cur_tp_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_tp_intron[bed_coords['start']]) <= entry_end
                        
            promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
            while promoterbool:
                cur_promoters.append(cur_promoter)
                cur_promoter = promoters.readline().strip().split('\t')
                promoterbool = len(cur_promoter) > 1 and cur_promoter[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_promoter[bed_coords['start']]) <= entry_end
                
            exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end
            while exonbool:
                cur_exons.append(cur_exon)
                cur_exon = exons.readline().strip().split('\t')
                exonbool = len(cur_exon) > 1 and cur_exon[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_exon[bed_coords['start']]) <= entry_end 

            intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 
            while intronbool:
                cur_introns.append(cur_intron)
                cur_intron = introns.readline().strip().split('\t')
                intronbool = len(cur_intron) > 1 and cur_intron[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_intron[bed_coords['start']]) <= entry_end 

            repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 
            while repeatbool:
                cur_repeats.append(cur_repeat)
                cur_repeat = repeats.readline().strip().split('\t')
                repeatbool = len(cur_repeat) > 1 and cur_repeat[bed_coords['chrom']]==entry_data[bed_coords['chrom']] and int(cur_repeat[bed_coords['start']]) <= entry_end 

            ## remove entries that are before the current entry here:
            for i in range(len(cur_fp_exons)-1, -1, -1):
                if cur_fp_exons[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_fp_exons[i][bed_coords['end']]) < entry_start:
                    cur_fp_exons.pop(i)
                    
            for i in range(len(cur_fp_introns)-1, -1, -1):
                if cur_fp_introns[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_fp_introns[i][bed_coords['end']]) < entry_start:
                    cur_fp_introns.pop(i)

            for i in range(len(cur_tp_exons)-1, -1, -1):
                if cur_tp_exons[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_tp_exons[i][bed_coords['end']]) < entry_start:
                    cur_tp_exons.pop(i)
                    
            for i in range(len(cur_tp_introns)-1, -1, -1):
                if cur_tp_introns[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_tp_introns[i][bed_coords['end']]) < entry_start:
                    cur_tp_introns.pop(i)
            
            for i in range(len(cur_promoters)-1, -1, -1):
                if cur_promoters[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_promoters[i][bed_coords['end']]) < entry_start:
                    cur_promoters.pop(i)

            for i in range(len(cur_exons)-1, -1, -1):
                if cur_exons[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_exons[i][bed_coords['end']]) < entry_start:
                    cur_exons.pop(i)

            for i in range(len(cur_introns)-1, -1, -1):
                if cur_introns[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_introns[i][bed_coords['end']]) < entry_start:
                    cur_introns.pop(i)

            for i in range(len(cur_repeats)-1, -1, -1):
                if cur_repeats[i][bed_coords['chrom']] != entry_data[bed_coords['chrom']] or int(cur_repeats[i][bed_coords['end']]) < entry_start:
                    cur_repeats.pop(i)                                        
                                    
            ## now compare the current entry to all these peaks
            ## here i am relying on the correct merging of the bed entries from the element files
            ## because i dont ever re-merge anything
            this_fp_exon_bp = 0.0
            this_fp_intron_bp = 0.0
            this_tp_exon_bp = 0.0
            this_tp_intron_bp = 0.0            
            this_promoter_bp = 0.0
            this_exon_bp = 0.0
            this_intron_bp = 0.0
            this_repeat_bp = 0.0
            
            for cfe in cur_fp_exons:
                this_fp_exon_bp += max(0, min(int(cfe[bed_coords['end']]), entry_end) - max(int(cfe[bed_coords['start']]), entry_start))
            for cfi in cur_fp_introns:
                this_fp_intron_bp += max(0, min(int(cfi[bed_coords['end']]), entry_end) - max(int(cfi[bed_coords['start']]), entry_start))
            for cte in cur_tp_exons:
                this_tp_exon_bp += max(0, min(int(cte[bed_coords['end']]), entry_end) - max(int(cte[bed_coords['start']]), entry_start))
            for cti in cur_tp_introns:
                this_tp_intron_bp += max(0, min(int(cti[bed_coords['end']]), entry_end) - max(int(cti[bed_coords['start']]), entry_start))
            for cp in cur_promoters:
                this_promoter_bp += max(0, min(int(cp[bed_coords['end']]), entry_end) - max(int(cp[bed_coords['start']]), entry_start))
            for ce in cur_exons:
                this_exon_bp += max(0, min(int(ce[bed_coords['end']]), entry_end) - max(int(ce[bed_coords['start']]), entry_start))
            for ci in cur_introns:
                this_intron_bp += max(0, min(int(ci[bed_coords['end']]), entry_end) - max(int(ci[bed_coords['start']]), entry_start))
            for cr in cur_repeats:
                this_repeat_bp += max(0, min(int(cr[bed_coords['end']]), entry_end) - max(int(cr[bed_coords['start']]), entry_start))
            ## add these overlaps to the summary trackers
            this_chr_fp_exon_bp += this_fp_exon_bp
            this_chr_fp_intron_bp += this_fp_intron_bp
            this_chr_tp_exon_bp += this_tp_exon_bp
            this_chr_tp_intron_bp += this_tp_intron_bp                                    
            this_chr_promoter_bp += this_promoter_bp
            this_chr_exon_bp += this_exon_bp
            this_chr_intron_bp += this_intron_bp
            this_chr_repeat_bp += this_repeat_bp
            
            total_fp_exon_bp += this_fp_exon_bp
            total_fp_intron_bp += this_fp_intron_bp
            total_tp_exon_bp += this_tp_exon_bp
            total_tp_intron_bp += this_tp_intron_bp            
            total_promoter_bp += this_promoter_bp
            total_exon_bp += this_exon_bp
            total_intron_bp += this_intron_bp
            total_repeat_bp += this_repeat_bp
            
            entry_out.write('\t'.join([entry_data[bed_coords['chrom']], entry_data[bed_coords['start']], entry_data[bed_coords['end']], str(this_fp_exon_bp), str(this_fp_exon_bp / this_length), str(this_fp_intron_bp), str(this_fp_intron_bp / this_length), str(this_tp_exon_bp), str(this_tp_exon_bp / this_length), str(this_tp_intron_bp), str(this_tp_intron_bp / this_length), str(this_promoter_bp), str(this_promoter_bp / this_length), str(this_exon_bp), str(this_exon_bp / this_length), str(this_intron_bp), str(this_intron_bp / this_length), str(this_repeat_bp), str(this_repeat_bp / this_length)])+'\n')

        summary_out.write('\t'.join(['genomewide', str(total_fp_exon_bp), str(total_fp_exon_bp / total_entry_bp), str(total_fp_intron_bp), str(total_fp_intron_bp / total_entry_bp), str(total_tp_exon_bp), str(total_tp_exon_bp / total_entry_bp), str(total_tp_intron_bp), str(total_tp_intron_bp / total_entry_bp), str(total_promoter_bp), str(total_promoter_bp/total_entry_bp), str(total_exon_bp), str(total_exon_bp/total_entry_bp), str(total_intron_bp), str(total_intron_bp/total_entry_bp), str(total_repeat_bp), str(total_repeat_bp/total_entry_bp)])+'\n')

    input_beds.close()
        
    end = time.clock()
    length = end - start
    print "Analysis complete, time: ", length        
        
if __name__=="__main__":
    # create the argument parser
    parser = argparse.ArgumentParser(description="Compute coverage statistics for each entry of a bed file. All input files should be sorted according to chromosome, strand (if applicable) and start position.")
    parser.add_argument("--fp_utr", help="The optional .bed file containing the full 5' UTR loci", default=None)
    parser.add_argument("--tp_utr", help="The optional .bed file containing the full 3' UTR loci", default=None)
    parser.add_argument("--full_utrs", nargs=4, help="4 .bed files containing the 5' UTR exons, 5' UTR introns, 3' UTR exons, and 3' UTR introns, in that order.", default=None, metavar=('FP_EXONS', 'FP_INTRONS', 'TP_EXONS', 'TP_INTRONS'))
    parser.add_argument("promoter_bed", help="The .bed file containing promoter loci")
    parser.add_argument("exon_bed", help="The .bed file containing exons")
    parser.add_argument("intron_bed", help="The .bed file containing introns")
    parser.add_argument("repeat_bed", help="The .bed file containing the repeat regions")
    parser.add_argument("input_bed", help="The input .bed file that you want to analyze")
    parser.add_argument("entrywise_output", help="The path to the desired output file containing each entry of the input bed file with the amount of overlap of each element")
    parser.add_argument("summary_output", help="The path to the desired output file containing summary statistics of the coverage of each type of element")
    pargs = parser.parse_args()
    
    if pargs.fp_utr and pargs.tp_utr:
        if (pargs.fp_utr and not pargs.tp_utr) or (pargs.tp_utr and not pargs.fp_utr):
            print "Need both 5' and 3' UTR files for UTR analysis; performing promoter/exon/intron/repeat/intergenic coverage only"
        else:
            compute_full_utr_coverage(pargs.fp_utr, pargs.tp_utr, pargs.promoter_bed, pargs.exon_bed, pargs.intron_bed, pargs.repeat_bed, pargs.input_bed, pargs.entrywise_output, pargs.summary_output)
    elif pargs.full_utrs:
        fp_exons, fp_introns, tp_exons, tp_introns = pargs.full_utrs
        compute_split_utr_coverage(fp_exons, fp_introns, tp_exons, tp_introns, pargs.promoter_bed, pargs.exon_bed, pargs.intron_bed, pargs.repeat_bed, pargs.input_bed, pargs.entrywise_output, pargs.summary_output)
    else:
        compute_coverage(pargs.promoter_bed, pargs.exon_bed, pargs.intron_bed, pargs.repeat_bed, pargs.input_bed, pargs.entrywise_output, pargs.summary_output)
