#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
import copy
import operator
#from scipy.interpolate import spline #this one is obsolete
from scipy.interpolate import BSpline #for interpolating coverage function
import matplotlib.patheffects


import sys
import itertools


#from CIAlign import FastaToArray


# python3 consensusSeq.py --infile /Users/lotti/Documents/Test_Case/gap_test/SevenJEV.afa


def FastaToArray(infile):
    '''
    Convert an alignment into a numpy array.
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

    nams = []
    seqs = []
    nam = ""
    seq = ""
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if line[0] == ">":
                    seqs.append(seq)
                    nams.append(nam)
                    seq = []
                    nam = line.replace(">", "")
            else:
                seq += list(line)
    seqs.append(seq)
    nams.append(nam)
    arr = np.array(seqs[1:])
    return (arr, nams[1:])

    
class Scale(matplotlib.patheffects.RendererBase):
    #Credits: Markus Piotrowski See: https://github.com/biopython/biopython/issues/850#issuecomment-225708297
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine=affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def findConsensus(alignment, consensus_type="majority"):
    '''
    '''
    consensus = []
    coverage = []
    numberOfSequences = len(alignment[:,0])

    #need the reverse of array to access every column
    for i in range(0,len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i], return_counts=True)
        count = dict(zip(unique, counts))
        unique_ng = unique[unique != "-"]
        counts_ng = counts[unique != "-"]
        count_ng = dict(zip(unique_ng, counts_ng))
        if '-' in count:
            nonGapContent = 1-(count['-']/numberOfSequences)
        else:
            nonGapContent = 1
        maxChar, maxCount = max(count.items(), key=operator.itemgetter(1))
        maxChar_ng, maxCount_ng = max(count_ng.items(), key=operator.itemgetter(1))

        # if there are an equal number of gap and non-gap characters at the
        # site, keep the non-gap character

        if maxCount_ng == maxCount or consensus_type == "majority_nongap":
            maxChar = maxChar_ng
        consensus.append(maxChar)
        coverage.append(nonGapContent)

    return consensus, coverage


def makePlot(consensus, coverage):

    x = np.arange(0, len(coverage), 1);
    y = coverage

    f = plt.figure()
    a = f.add_subplot('311')
    a.plot(x,y)
    a.get_xaxis().set_visible(False)

    b = f.add_subplot('312')

    xnew = np.linspace(x.min(),x.max(),300) #300 represents number of points to make between T.min and T.max

    power_smooth = spline(x,y,xnew)

    b.plot(xnew,power_smooth)
    b.get_xaxis().set_visible(False)
    f.savefig('/Users/lotti/blub.png')

    #c = f.add_subplot('313')
    #txt = plt.text(x,consensus)
    #txt.set_path_effects([Scale(x_scale, y_scale)])
    plt.figure()
    txt = plt.text(0, 0, "T", fontsize=64,color='red')
    txt.set_path_effects([Scale(1,5)])
    txt = plt.text(0.1, 0, "G", fontsize=64,color='green')
    txt.set_path_effects([Scale(1,3)])
    txt = plt.text(0.2, 0, "B", fontsize=64,color='blue')
    txt.set_path_effects([Scale(1,7)])
    plt.show()

def main():
    # this is just for testing purposes
    print('blaaa')
    parser = argparse.ArgumentParser(
            description='''Improve a multiple sequence alignment''')

    parser.add_argument("--infile", dest='infile', type=str,
                        help='path to input alignment')

    args = parser.parse_args()

    arr, orig_nams = FastaToArray(args.infile)

    consensus, coverage = findConsensus(arr)
    makePlot(consensus, coverage)
    print(consensus)
    print(coverage)



if __name__ == "__main__":
    main()
