#! /usr/bin/env python

import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
import copy
import operator

from CIAlign import FastaToArray
from consensusSeq import findConsensus

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
