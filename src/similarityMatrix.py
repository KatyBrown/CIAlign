#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import itertools

def calculateSimilarityMatrix(arr, nams, minoverlap=1,
                              keepgaps=False, outfile=None, dp=4):
    '''
    Calculates a pairwise similarity matrix for the alignment
    dp = decimal places (in the output file only)
    outfile = path to outfile
    keepgaps = should positions with gaps in either sequence in the pair
    be discarded before calculating identity?
    minoverlap = minimum number of positions at which both sequences should
    have a non-gap character for the percentage to not be 0
    '''
    ident = np.empty((len(arr), len(arr)))
    for i, j in itertools.combinations_with_replacement(range(len(arr)), 2):
        p1 = arr[i]
        p2 = arr[j]
        if not keepgaps:
            nongap = (p1 != "-") & (p2 != "-")
        else:
            nongap = np.array([True] * len(p1))
        matches = sum(p1[nongap] == p2[nongap])
        length = sum(nongap)
        # minimum overlap is currently used regardless of the keepgaps settings
        if length >= minoverlap and sum((p1 != "-") & (p2 != "-")) >= minoverlap:
            perc = matches / length
        else:
            perc = 0
        ident[i, j] = perc
        ident[j, i] = perc

    if outfile:
        out = open(outfile, "w")
        out.write("\t%s\n" % ("\t".join(nams)))
        for i, line in enumerate(ident):
            out.write("%s\t%s\n" % (
                    nams[i], "\t".join([str(round(L, dp)) for L in line])))
        out.close()
    return(ident)