#! /usr/bin/env python

import logging
import argparse
import numpy as np

def FastaToDict(infile):
    '''
    Converts a fasta file to a dictionary
    
    Keys are the sequence names (without ">") and values are the
    sequences.
    
    Parameters
    ----------
    infile: str
        path to fasta file
        
    Returns
    -------
    dict
        dictionary where keys are sequence names and values are sequence or
        sequence lengths
    '''
    D = dict()
    seq = []
    nam = ""
    nams = []
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) != 0:
                if line[0] == ">":
                    if len(seq) != 0:
                        seq = "".join(seq)
                        D[nam] = seq
                        nams.append(nam)
                        seq = []
                    nam = line.replace(">", "")
                else:
                    seq.append(line)
    seq = "".join(seq)
    D[nam] = seq
    nams.append(nam)
    return (D, nams)


def alignToArray(fasta_dict):
    '''
    Convert an alignment dictionary into a numpy array.
    Parameters
    ----------
    fasta_dict: dict
        dictionary based on a fasta file with sequence names as keys and
        sequences as values
        
    Returns
    -------
    arr: np.array
        2D numpy array in the same order as fasta_dict where each row
        represents a single column in the alignment and each column a
        single sequence.
    '''
    arr = []
    for nam in sorted(fasta_dict.keys()):
        arr.append(np.array(list(fasta_dict[nam])))
    arr = np.array(arr)
    return (arr)


def removeIndels(arr, log, sliding_window_size=10, mincov=10):
    log.info("Removing indels\n")
    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column 
    sums = sum(boolarr)

    # run a sliding window along the alignment and check for regions
    # which have higher coverage at the ends of the window than in the
    # middle - store these in put_indels
    put_indels = set()
    for size in range(9, sliding_window_size):
        for i in range(0, len(sums)+1 - size):
            these_sums = sums[i:i+size]
            # take the number of non gap positions
            # for each column in this window
            ns = np.array(range(i, i+size))
            left = these_sums[0]
            right = these_sums[-1]
            # record sites with lower coverage than their flanking seqs
            # and less than mincov total coverage
            x = set(ns[(these_sums < left) &
                       (these_sums < right) &
                       (these_sums < mincov)])
            put_indels = put_indels | x
    put_indels = np.array(sorted(list(put_indels)))
    # for the putative indels, check if any sequence without the indel
    # flanks the indel site - if it does it confirms it is an indel
    rmpos = []
    for p in put_indels:
        has_flanks = 0
        for a in arr:
            thispos = a[p]
            if thispos == "-":
                leftsum = sum(a[:p] != "-")
                rightsum = sum(a[p:] != "-")
                if leftsum > mincov and rightsum > mincov:
                    has_flanks += 1
        if has_flanks != 0:
            rmpos.append(p)
    # make a list of positions to keep
    keeppos = np.arange(0, len(sums))
    keeppos = np.invert(np.in1d(keeppos, rmpos))
    arr = arr[:, keeppos]


def writeOutfile(outfile, arr, fasta_dict, orig_nams):
    out = open(outfile, "w")
    F = dict()
    for i, nam in enumerate(sorted(fasta_dict.keys())):
        F[nam] = "".join(list(arr[i]))
        
    # ensures the output file is in the same order as the input file
    for nam in orig_nams:
        out.write(">%s\n%s\n" % (nam, F[nam]))
    out.close()


def main():
    parser = argparse.ArgumentParser(
            description='''Improve a multiple sequence alignment''')

    parser.add_argument("--infile", dest='infile', type=str,
                        help='path to input alignment')
    parser.add_argument("--outfile", dest='outfile', type=str, default="",
                        help="path to output alignment")
    parser.add_argument("--sliding_window_size_indels", dest="slidingwindowsize",
                        type=int, default=200,
                        help="sliding window size to remove indels")
    parser.add_argument("--max_coverage_indels", dest="maxcov",
                        type=int, default=20,
                        help="maximum coverage of removed indels")
    parser.add_argument("--remove_indels", dest="rmi",
                        action="store_true")

    args = parser.parse_args()

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    
    handler = logging.FileHandler("log.log")
    handler.setLevel(logging.INFO)
    
    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    
    # add the handlers to the logger
    log.addHandler(handler)
    
    fasta_dict, orig_nams = FastaToDict(args.infile)
    arr = alignToArray(fasta_dict)

    if args.rmi:
        removeIndels(arr, log, args.slidingwindowsize, args.maxcov)
    
    writeOutfile(args.outfile, arr, fasta_dict, orig_nams)


    

if __name__ == "__main__":
    main()