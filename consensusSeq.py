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
import os
from matplotlib.font_manager import FontProperties
#from scipy.interpolate import spline #this one is obsolete
import matplotlib.patheffects

#from PIL import ImageFont


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
    return {'A': '#1ed30f',
            'G': '#f4d931',
            'T': '#f43131',
            'C': '#315af4',
            'N': '#b2b2b2',
            "-": '#FFFFFF',
            "U": '#f43131'}


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
    return (rect_height_points * 1.8)


def PixelsToPoints(pixels, dpi):
    return (pixels * (72 / dpi))


def getLetters(typ='nt', fontname='monospace', dpi=500):
    if typ == 'nt':
        colours = getNtColours()
    elif typ == 'aa':
        colours = getAAColours()
    for base in colours.keys():
        f = plt.figure(figsize=(1, 1), dpi=dpi, edgecolor='black')
        a = f.add_subplot(111)
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        fs = getFontSize(f, a, 1)
        a.text(0.5, 0, base, fontsize=fs, fontdict={'family': 'monospace',
                                                  'name': fontname},
               color=colours[base], va='baseline', ha='center')
        plt.gca().set_axis_off()
        a.margins(0, 0)
        f.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, wspace=None, hspace=None)
        a.set_frame_on(False)
        f.savefig("%s_temp.png" % base, dpi=500,
                  pad_inches=0)

# class Scale(matplotlib.patheffects.RendererBase):
#     #Credits: Markus Piotrowski See: https://github.com/biopython/biopython/issues/850#issuecomment-225708297
#     def __init__(self, sx, sy=None):
#         self._sx = sx
#         self._sy = sy
#
#     def draw_path(self, renderer, gc, tpath, affine, rgbFace):
#         affine=affine.identity().scale(self._sx, self._sy)+affine
#         renderer.draw_path(gc, tpath, affine, rgbFace)


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


def sequence_logo(alignment,
                  figname,
                  typ='nt',
                  figfontname='Arial',
                  figdpi=500):
    for seq in alignment:
        if np.any(seq == "U"):
            that_letter = "U"
        if np.any(seq == "T"):
            that_letter = "T"
    f = plt.figure(figsize=(10, 2), dpi=figdpi)
    a = f.add_subplot(111)
    getLetters(typ=typ, fontname=figfontname, dpi=figdpi)
    a.set_xlim(0, len(alignment[0,:]))
    a.set_ylim(0, 3.1)
    limits = a.axis()
    getLetters(typ=typ, fontname=figfontname)
    for i in range(0, len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i],
                                   return_counts=True)
        count = dict(zip(unique, counts))
        height_per_base, info_per_base = calc_entropy(count,
                                                      len(alignment[:,0]),
                                                      that_letter=that_letter)

        height_sum_higher = 0
        for base, height in height_per_base.items():
            if height > 0:
                L = plt.imread("%s_temp.png" % base)
                print (base, height)
                a.imshow(L, extent=(i, i+1, height_sum_higher, height_sum_higher+height))
                height_sum_higher += height
    a.axis(limits)
    f.set_size_inches(10, 2)
    #a.set_xticks(np.arange(0, len(alignment[0,:])))
    #a.set_xticklabels(np.arange(1, len(alignment[0,:]) + 1))
    #a.set_yticks(np.arange(0, 3.1, 1))
    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.set_xlabel("Position")
    a.set_ylabel("Bit Score")
    if typ == 'nt':
        allbases = getNtColours()
    elif typ == 'aa':
        allbases = getAAColours()
    for base in allbases:
        os.unlink("%s_temp.png" % base)
    f.savefig(figname, dpi=figdpi, bbox_inches='tight')

"""
    #plt.figure(figsize=(len(alignment[0,:]),2.5))

    #plt.xkcd()
    axes = plt.gca()
    axes.set_xlim([0,len(alignment[0,:])+1])
    axes.set_ylim([0,3])
    seq_count = len(alignment[:,0])
    x = 0
    font = FontProperties()
    font.set_size(80)
    #print("fuck die henne", font.ascent())

    for i in range(0,len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i], return_counts=True)
        count = dict(zip(unique, counts))
        height_per_base, info_per_base = calc_entropy(count, seq_count)
        print('height', height_per_base)

        height_sum_higher = 0
        for base, height in height_per_base.items():
            if height > 0:
                txt = plt.text(x, height_sum_higher, base, fontsize=70, color='red')
                #txt = plt.text(x, height_sum_higher, base, fontsize=80, color='red')
                #txt.set_path_effects([Scale(1,height)])
                height_sum_higher += height
        # coordinates and then scale aha aha aha
        # txt = plt.text(0, 0, "T", fontsize=64,color='red')
        # txt.set_path_effects([Scale(1,5)])
        # txt = plt.text(0, 0.1, "G", fontsize=64,color='green')
        # txt.set_path_effects([Scale(1,3)])
        # txt = plt.text(0.2, 0, "B", fontsize=64,color='blue')
        # txt.set_path_effects([Scale(1,7)])
        x += 1

    plt.show()
    plt.savefig('plotileini.png')
"""

def sequence_bar_logo(alignment,
                      figname,
                      typ='nt',
                      figdpi=500):

    for seq in alignment:
        print(type(seq))
        if np.any(seq == "U"):
            that_letter = "U"
        if np.any(seq == "T"):
            that_letter = "T"
    # if len(alignment[0,:]) < 65536:
    #     plt.figure(figsize=(len(alignment[0,:]) + 1,4))
    # else:
    plt.figure(1, figsize=(len(alignment[0,:]), 4), frameon=False,
               dpi=figdpi)
    fig, ax = plt.subplots()

    #plt.xkcd()
    axes = plt.gca()
    axes.set_xlim([-0.5, len(alignment[0,:])-0.5])
    axes.set_ylim([0, 3.1])
    seq_count = len(alignment[:,0])
    x = 0
    width = 0.75
    #todo alternative color scheme
    ind = np.arange(len(alignment[0,:]))
    A_height = []
    G_height = []
    C_height = []
    U_height = []

    for i in range(0,len(alignment[0,:])):
        unique, counts = np.unique(alignment[:,i], return_counts=True)
        count = dict(zip(unique, counts))
        height_per_base, info_per_base = calc_entropy(count, seq_count, that_letter)
        print('height', height_per_base)

        A_height.append(height_per_base["A"])
        G_height.append(height_per_base["G"])
        C_height.append(height_per_base["C"])
        U_height.append(height_per_base[that_letter])

    if typ == 'nt':
        colours = getNtColours()
    elif typ == 'aa':
        colours = getAAColours()

    bar_plot = plt.bar(ind, A_height, width, color=colours['A'])
    # for idx,rect in enumerate(bar_plot):
    #     height = rect.get_height()
    #     ax.text(rect.get_x() + rect.get_width(), height, "A", ha='center', va='bottom')

    plt.bar(ind, G_height, width, bottom=A_height, color=colours['G'])
    plt.bar(ind, C_height, bottom=[i+j for i,j in zip(A_height, G_height)], width=width, color=colours['C'])
    plt.bar(ind, U_height, bottom=[i+j+k for i,j,k in zip(A_height, G_height, C_height)], width=width,
                                   color=colours['U'])
    #plt.xticks(np.arange(0, len(alignment[0,:])), np.arange(1, len(alignment[0,:]) + 1))
    #plt.yticks(np.arange(0, 3.1, 1))
    plt.xlabel("Position")
    plt.ylabel("Bit Score")

    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)

    plt.savefig(figname, bbox_inches='tight', dpi=figdpi)
    #todo tidy up and use normal names and normals files



def calc_entropy(count, seq_count, that_letter):

    # total number of Sequences - gap number
    # adjust total height later to make up for gaps - i think that's covered (?)

    sample_size_correction = (1/log(2)) * (3/(2*seq_count)) # wahh here or later?
    gap_correction = seq_count
    if count.get("-"):
        seq_count -= count.get("-")
    gap_correction = seq_count/gap_correction
    print('gap correction', gap_correction)
    info_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    freq_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    height_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    entropy_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    entropy = 0
    if seq_count == 0:
        return height_per_base, info_per_base

    for base, quantity in count.items():
        if base != "-":
            frequency = quantity/seq_count
            freq_per_base[base] = frequency
            entropy -= frequency*log(frequency, 2)
            info_per_base[base] = 2 + frequency*log(frequency, 2)
            entropy_per_base[base] = -frequency*log(frequency,2)
    information_per_column = 2-entropy-sample_size_correction
    print("info", info_per_base)
    for base, quantity in info_per_base.items():
        if freq_per_base[base]*information_per_column < 0:
            height_per_base[base] = 0
        else:
            #scale to accomodate gaps, or do we even need to that here if already
            #used in the sample size correction?
            height_per_base[base] = gap_correction*freq_per_base[base]*information_per_column

    return height_per_base, info_per_base

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
    sequence_bar_logo(arr)
    makePlot(consensus, coverage)
    print(consensus)
    print(coverage)



if __name__ == "__main__":
    main()
