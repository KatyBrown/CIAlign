#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from math import log
import logging
import argparse
import numpy as np
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import copy
import operator
#from scipy.interpolate import spline #this one is obsolete
import matplotlib.patheffects


import sys
import itertools


#from CIAlign import FastaToArray


# python3 consensus Seq.py --infile /Users/lotti/Documents/Test_Case/gap_test/SevenJEV.afa


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


def getAAColours():
    return {'D':'#E60A0A',
            'E':'#E60A0A',
            'C':'#E6E600',
            'M':'#E6E600',
            'K':'#145AFF',
            'R':'#145AFF',
            'S':'#FA9600',
            'T':'#FA9600',
            'F':'#3232AA',
            'Y':'#3232AA',
            'N':'#00DCDC',
            'Q':'#00DCDC',
            'G':'#EBEBEB',
            'L':'#0F820F',
            'V':'#0F820F',
            'I':'#0F820F',
            'A':'#C8C8C8',
            'W':'#B45AB4',
            'H':'#8282D2',
            'P':'#DC9682',
            'X': '#b2b2b2',
            '-': '#FFFFFF00'}


def getNtColours():
    return {'A': '#f43131',
            'G': '#f4d931',
            'T': '#315af4',
            'C': '#1ed30f',
            'N': '#b2b2b2',
            "-": '#FFFFFF',
            "U": '#315af4'}


def getAxisUnits(figure, subplot):
    axis_dimensions = subplot.transData.transform([(subplot.get_xlim()[1], subplot.get_ylim()[1]),(0, 0)])- subplot.transData.transform((0,0))
    D = dict()
    D['axis_height_px'] = axis_dimensions[0][1]
    D['axis_width_px'] = axis_dimensions[0][0]
    D['axis_bottom'], D['axis_top'] = subplot.get_ylim()
    D['axis_left'], D['axis_right'] = subplot.get_xlim()
    D['axis_width_u'] = D['axis_right'] - D['axis_left']
    D['axis_height_u'] = D['axis_top'] - D['axis_bottom']
    D['u_height_px'] = D['axis_height_px'] / D['axis_height_u']
    D['u_width_px'] =  D['axis_width_px'] / D['axis_width_u']
    return (D)


def getFontSize(figure, subplot, rectangle_height_u):
    D = getAxisUnits(figure, subplot)
    rect_perc_height = rectangle_height_u / (D['axis_top'] - D['axis_bottom'])
    rect_height_pixels = D['axis_height_px'] * rect_perc_height
    dpi = figure.dpi
    rect_height_points = PixelsToPoints(rect_height_pixels, dpi)
    return (rect_height_points * 1.4)


def PixelsToPoints(pixels, dpi):
    return (pixels * (72 / dpi))


def getLetters(typ='nt', fontfamily='monospace'):
    if typ == 'nt':
        colours = getNtColours()
    elif typ == 'aa':
        colours = getAAColours()
    D = dict()
    for base in colours.keys():
        f = plt.figure(figsize=(1, 1), dpi=500)
        a = f.add_subplot(111)
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        fs = getFontSize(f, a, 1)
        a.text(0, 0, base, fontsize=fs, fontdict={'family': fontfamily},
               color=colours[base])
        a.set_axis_off()
        f.canvas.draw()
        letter = np.frombuffer(f.canvas.tostring_rgb(), dtype=np.uint8)
        letter = letter.reshape(f.canvas.get_width_height()[::-1] + (3,))
        D[base] = letter
    return (D)

def tempPlotLetters(string, heights,
                    typ='nt',
                    figheight=2,
                    figwidth=10,
                    figfontfamily='monospace'):
    f = plt.figure(figsize=(figwidth, figheight), dpi=500)
    a = f.add_subplot(111)
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    letters = getLetters(typ=typ)
    D = getAxisUnits(f, a)
    y_int = 1 / len(string)
    x = 0
    y = 0
    for char in string:
        a.matshow(letters[char], extent=(x, x+0.2, 0, y), resample=False)
        x += 0.2
        y += y_int    
    f.set_size_inches(figwidth, figheight)
    f.savefig("test.png", dpi=500)

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
        if counts_ng.size == 0:
            count_ng = {"N": len(alignment[:,i])}
            nonGapContent = 0
        else:
            count_ng = dict(zip(unique_ng, counts_ng))
            if '-' in count:
                nonGapContent = 1-(count['-']/numberOfSequences)
            else:
                nonGapContent = 1

        # dealing with gap only collumns
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
    print(len(x), len(y))

    f = plt.figure()
    a = f.add_subplot('311')
    a.plot(x,y)
    a.get_xaxis().set_visible(False)

    b = f.add_subplot('312')

    #xnew = np.linspace(x.min(),x.max(),300) #300 represents number of points to make between T.min and T.max

    t, c, k = interpolate.splrep(x, y, s=0, k=4)
    N = 100
    xmin, xmax = x.min(), x.max()
    xx = np.linspace(xmin, xmax, N)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    b.plot(xx, spline(xx))
    b.get_xaxis().set_visible(False)
    f.savefig('blub.png')

def sequence_logo(alignment):

    plt.figure()
    seq_count = len(alignment[:,0])

    for i in range(0,len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i], return_counts=True)
        count = dict(zip(unique, counts))
        entropy, info_per_base, freq_per_base, height_per_base = calc_entropy(count, seq_count)

        txt = plt.text(0, 0, "T", fontsize=64,color='red')
        txt.set_path_effects([Scale(1,5)])
        txt = plt.text(0, 0.1, "G", fontsize=64,color='green')
        txt.set_path_effects([Scale(1,3)])
        txt = plt.text(0.2, 0, "B", fontsize=64,color='blue')
        txt.set_path_effects([Scale(1,7)])

    plt.show()
    plt.savefig('plotileini.png')

def calc_entropy(count, seq_count):

    # total number of Sequences - gap number
    # adjust total height later to make up for gaps
    if count.get("-"):
        seq_count -= count.get("-")
    info_per_base = {"A": 0, "G": 0, "U": 0, "C": 0}
    freq_per_base = info_per_base
    height_per_base = info_per_base
    entropy = 0
    if seq_count == 0:
        return entropy, info_per_base, freq_per_base, height_per_base

    for base, quantity in count.items():
        if base != "-":
            frequency = quantity/seq_count
            freq_per_base[base] = frequency
            entropy -= frequency*log(frequency, 2)
            info_per_base[base] = 2 + frequency*log(frequency, 2)
    for base, quantity in info_per_base.items():
        height_per_base[base] = freq_per_base[base]*(2-entropy)
    print('freq', freq_per_base)
    print('height', height_per_base)

    return entropy, info_per_base, freq_per_base, height_per_base

def main():
    # this is just for testing purposes
    print('consensus test')
    parser = argparse.ArgumentParser(
            description='''Improve a multiple sequence alignment''')

    parser.add_argument("--infile", dest='infile', type=str,
                        help='path to input alignment')

    args = parser.parse_args()

    arr, orig_nams = FastaToArray(args.infile)

    consensus, coverage = findConsensus(arr)
    sequence_logo(arr)
    makePlot(consensus, coverage)
    print(consensus)
    print(coverage)



if __name__ == "__main__":
    main()
