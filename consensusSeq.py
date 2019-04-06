#! /usr/bin/env python

import logging
import argparse
import numpy as np
#import matplotlib.pyplot as plt
import copy
import operator

# python3 consensusSeq.py --infile /Users/lotti/Documents/Test_Case/gap_test/SevenJEV.afa

from CIAlign import FastaToDict, DictToArray

def findConsensus(alignment):

    consensus = []

    #need the reverse of array to access every column
    for i in range(0,len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i], return_counts=True)
        count = dict(zip(unique, counts))
        print(count)
        maxChar = max(count.items(), key=operator.itemgetter(1))[0]
        consensus.append(maxChar)

    return consensus

def main():
    print('blaaa')
    parser = argparse.ArgumentParser(
            description='''Improve a multiple sequence alignment''')

    parser.add_argument("--infile", dest='infile', type=str,
                        help='path to input alignment')

    args = parser.parse_args()

    fasta_dict, orig_nams = FastaToDict(args.infile)

    # convert the fasta file dictionary into a numpy array
    arr = DictToArray(fasta_dict)

    consensus = findConsensus(arr)
    print(len(consensus))
    print(consensus)



if __name__ == "__main__":
    main()
