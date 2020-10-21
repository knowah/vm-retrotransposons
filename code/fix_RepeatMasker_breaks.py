#!/usr/bin/env python3

import argparse
import sys
import gzip
import contextlib

RM_FIELDS = ['chrom', 'start', 'end', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'element_ID']
class RM_Entry(object):
    __slots__ = RM_FIELDS
    def __str__(self):
        return "\t".join([str(self.__getattribute__(x)) for x in RM_FIELDS])

def process_line(ln):
    tokens = ln.rstrip('\n').split('\t')
    # chrom start end strand repName repClass repFamily repStart repEnd element.ID
    entry = RM_Entry()
    entry.chrom  = tokens[5]
    entry.start  = int(tokens[6])
    entry.end    = int(tokens[7])
    entry.strand = tokens[9]
    entry.repName   = tokens[10]
    entry.repClass  = tokens[11]
    entry.repFamily = tokens[12]
    entry.repStart = int(tokens[13])
    entry.repEnd   = int(tokens[14])
    entry.element_ID = int(tokens[16])
    return entry

def process_line_alt(ln):
    tokens = ln.rstrip('\n').split('\t')
    #chrom chromStart chromEnd name score strand swScore milliDiv milliDel milliIns genoLeft repClass repFamily repStart repEnd repLeft
    entry = RM_Entry()
    entry.chrom  = tokens[0]
    entry.start  = int(tokens[1])
    entry.end    = int(tokens[2])
    entry.strand = tokens[5]
    entry.repName   = tokens[3]
    entry.repClass  = tokens[11]
    entry.repFamily = tokens[12]
    entry.repStart = abs(int(tokens[13]))
    entry.repEnd   = int(tokens[14])
    entry.element_ID = lineno
    return entry


@contextlib.contextmanager
def fopen(filename=None, mode='rt'):
    if filename and filename != "-":
        openfn = open if not filename.endswith('.gz') else gzip.open
        f = openfn(filename, mode)
    else: 
        f = sys.stdin if 'r' in mode else sys.stdout
    
    try:
        yield f
    finally:
        if f is not sys.stdin and f is not sys.stdout:
            f.close()

def fix_RM_breaks(infile, breaks=500000, seq_gap=0, rep_gap=1, alt_format=False, outfile=None):
    # process RepeatMasker file by line
    read_entry = process_line if not alt_format else process_line_alt
    
    with fopen(infile, 'rt') as inf, fopen(outfile, 'wt') as outf:
        lineno = 1

        firstln = inf.readline()
        while firstln.startswith('#'):
            # skip past header lines
            firstln = inf.readline()

        prev = read_entry(firstln)
        merge_elements = False
        elem_ids = (None, None) # when merge_elements == True, this stores the element IDs to merge
    
        for line in inf:
            lineno += 1
            curr = read_entry(line)
    
            # fail if input file not sorted by start
            if prev.chrom == curr.chrom and curr.start < prev.start:
                raise IOError("ERROR: Entries out of order (line {})".format(lineno))
    
            # check if this element and the previous one need to be merged
            if not merge_elements \
               and curr.start % breaks <= seq_gap \
               and curr.chrom == prev.chrom \
               and curr.strand == prev.strand \
               and curr.start - prev.end <= seq_gap \
               and curr.repName == prev.repName \
               and curr.repClass == prev.repClass \
               and curr.repFamily == prev.repFamily \
               and ((curr.strand == "+" and curr.repStart - prev.repEnd <= rep_gap) or \
                    (curr.strand == "-" and prev.repStart - curr.repEnd <= rep_gap)):
                merge_elements = True
                elem_ids = (prev.element_ID, curr.element_ID)
            
            # merge elements if the element ID matches the ID to be changed
            if merge_elements and curr.element_ID == elem_ids[1]:
                curr.element_ID = elem_ids[0]
            else:
                merge_elements = False
            
            outf.write(str(prev)+"\n") # output new entry for the previous line
            prev = curr
        
        outf.write(str(prev)+"\n") # print last line

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile",  type=str, help="Input (tab-separated) RepeatMasker file (must be sorted)")
    parser.add_argument("breaks",  type=int, default=500000, help="Breakpoints to merge upon [default: 500000]")
    parser.add_argument("-o", "--outfile", type=str, help="Output file (leave blank or - for stdout)")
    parser.add_argument("-s", "--seq_gap", type=int, default=0, help="Maximum (reference sequence) gap in bp between broken subelements")
    parser.add_argument("-r", "--rep_gap", type=int, default=1, help="Maximum gap between repEnd & repStart of broken subelements")
    parser.add_argument("--alt_format", action="store_true", help="Use alternative input format (e.g. from mouse strain genomes)")
    args = parser.parse_args()
    
    fix_RM_breaks(args.infile, args.breaks, args.seq_gap, args.rep_gap, args.alt_format, args.outfile)


