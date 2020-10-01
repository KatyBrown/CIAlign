#!/usr/bin/env python3
import numpy as np
import matplotlib
import itertools
matplotlib.use('Agg')


def calculateSimilarityMatrix(arr, nams, minoverlap=1,
                              keepgaps=0, outfile=None, dp=4):
    '''
    Calculates a pairwise similarity matrix for the alignment - for each pair
    of sequences, number of matching bases / total number of bases.
    Outputs a text file containing this matrix.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    minoverlap: int
        Minimum overlap between two sequences to have non-zero similarity
        in the similarity matrix.
    keepgaps: int
        Include positions with gaps in the
        similarity matrix calculation. Can be 0 - exclude positions which
        are gaps in either or both sequences, 1 - exclude positions which are
        gaps in both sequences, 2 - consider all positions regardless of gaps
    outfile: str
        Path to output file.  Can be None (for no output)
    dp: int
        Number of decimal places for the similarity matrix (output file only).

    Returns
    -------
    ident: np.array
        The identity matrix as a numpy array
    '''
    ident = np.empty((len(arr), len(arr)))
    for i, j in itertools.combinations_with_replacement(range(len(arr)), 2):
        p1 = arr[i]
        p2 = arr[j]
        if keepgaps == 0:
            # exclude positions which are gaps in either or both seqs
            nongap = (p1 != "-") & (p2 != "-")
        elif keepgaps == 1:
            # exclude positions which are gaps in both seqs
            nongap = (p1 != "-") | (p2 != "-")
        elif keepgaps == 2:
            # don't exclude anything
            nongap = np.array([True] * len(p1))
        matches = sum(p1[nongap] == p2[nongap])
        length = sum(nongap)
        # minimum overlap is currently used regardless of the keepgaps setting
        if length >= minoverlap and sum(
                (p1 != "-") & (p2 != "-")) >= minoverlap:
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
