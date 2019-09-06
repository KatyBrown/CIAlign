#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from math import log
import logging
import argparse
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import operator
import os
#from scipy.interpolate import spline #this one is obsolete
import matplotlib.patheffects
import math
from matplotlib import gridspec
import CIAlign.utilityFunctions as utilityFunctions

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
        colours = utilityFunctions.getNtColours()
    elif typ == 'aa':
        colours = utilityFunctions.getAAColours()
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
        f.savefig("%s_temp.tiff" % base, dpi=500,
                  pad_inches=0)
        plt.close()

# class Scale(matplotlib.patheffects.RendererBase):
#     #Credits: Markus Piotrowski See: https://github.com/biopython/biopython/issues/850#issuecomment-225708297
#     def __init__(self, sx, sy=None):
#         self._sx = sx
#         self._sy = sy
#
#     def draw_path(self, renderer, gc, tpath, affine, rgbFace):
#         affine=affine.identity().scale(self._sx, self._sy)+affine
#         renderer.draw_path(gc, tpath, affine, rgbFace)


def findConsensus(alignment, log, consensus_type="majority"):
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


def makeCoveragePlot(coverage, dest, dpi=300, height=3, width=5,
                     colour='#007bf5'):
    fontsize = 1500 / dpi
    x = np.arange(0, len(coverage), 1);
    y = coverage
    xmin, xmax = x.min(), x.max()
    N = 100

    xx = np.linspace(xmin, xmax, N)

    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot('211')
    a.plot(x, y, color=colour)
    a.set_xlabel('Position', fontsize=fontsize)
    a.set_ylabel('Coverage (Raw)', fontsize=fontsize)
    a.set_xticks([0, xmax])
    a.set_xticklabels([0, xmax], fontsize=fontsize)
    a.set_yticks(np.arange(0, 1.1, 0.5))
    a.set_yticklabels(np.arange(0, 1.1, 0.5), fontsize=fontsize)
    b = f.add_subplot('212')

    # polynomial interpolation
    #c = f.add_subplot('313')
    #z = np.polyfit(x, bla, 30)
    #p = np.poly1d(z)
    #c.plot(xx, p(xx))
    #xnew = np.linspace(x.min(),x.max(),300) #300 represents number of points to make between T.min and T.max

    t, c, k = interpolate.splrep(x, y, s=0, k=4)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    b.plot(xx, spline(xx), color=colour)
    b.set_xlabel('Position', fontsize=fontsize)
    b.set_ylabel('Coverage (Smoothed)', fontsize=fontsize)
    b.set_xticks([0, xmax])
    b.set_xticklabels([0, xmax], fontsize=fontsize)
    b.set_yticks(np.arange(0, 1.1, 0.5))
    b.set_yticklabels(np.arange(0, 1.1, 0.5), fontsize=fontsize)
    f.savefig(dest, dpi=dpi, bbox_inches='tight')


def sequence_logo(alignment,
                  figname,
                  typ='nt',
                  figfontname='Arial',
                  figdpi=300,
                  figrowlength=50):
    for seq in alignment:
        if np.any(seq == "U"):
            that_letter = "U"
        if np.any(seq == "T"):
            that_letter = "T"
    alignment_width = len(alignment[0,:])
    if alignment_width < figrowlength:
        figrowlength = alignment_width
    nsegs = math.ceil(alignment_width / figrowlength)
    f = plt.figure(figsize=(figrowlength, nsegs*2), dpi=figdpi)
    gs = gridspec.GridSpec(ncols=1, nrows=nsegs)
    getLetters(typ=typ, fontname=figfontname, dpi=figdpi)
    rstart = 0
    rend = rstart + figrowlength
    for n in range(nsegs):
        if rend > alignment_width:
            rend = alignment_width
        a = plt.subplot(gs[n])
        a.set_xlim(rstart, rstart+figrowlength)
        a.set_ylim(0, 3.1)
        limits = a.axis()

        for i in range(rstart, rend):

            unique, counts = np.unique(alignment[:,i],
                                       return_counts=True)
            count = dict(zip(unique, counts))
            height_per_base, info_per_base = calc_entropy(count,
                                                          len(alignment[:,0]),
                                                          that_letter=that_letter,
                                                          typ=typ)

            height_sum_higher = 0
            for base, height in height_per_base.items():
                if height > 0:
                    L = plt.imread("%s_temp.tiff" % base)
                    a.imshow(L, extent=(i, i+1, height_sum_higher, height_sum_higher+height),
                             filternorm=False)
                    height_sum_higher += height
        a.axis(limits)
        a.set_xticks([rstart, rend])
        a.set_xticklabels([rstart, rend])
        a.set_yticks(np.arange(0, 3.1, 1))
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        if n == (nsegs - 1):
            a.set_xlabel("Position")
        a.set_ylabel("Bit Score")
        rstart += figrowlength
        rend += figrowlength

    if typ == 'nt':
        allbases = utilityFunctions.getNtColours()
    elif typ == 'aa':
        allbases = utilityFunctions.getAAColours()
    for base in allbases:
        os.unlink("%s_temp.tiff" % base)
    f.savefig(figname, dpi=figdpi, bbox_inches='tight')
    plt.close()


def sequence_bar_logo(alignment,
                      figname,
                      typ='nt',
                      figdpi=500,
                      figrowlength=50):

    for seq in alignment:
        # print(type(seq))
        if np.any(seq == "U"):
            that_letter = "U"
        if np.any(seq == "T"):
            that_letter = "T"
    # if len(alignment[0,:]) < 65536:
    #     plt.figure(figsize=(len(alignment[0,:]) + 1,4))
    # else:

    alignment_width = len(alignment[0,:])
    if alignment_width < figrowlength:
        figrowlength = alignment_width
    nsegs = math.ceil(alignment_width / figrowlength)
    f = plt.figure(figsize=(figrowlength/5, nsegs*2), dpi=figdpi)
    gs = gridspec.GridSpec(ncols=1, nrows=nsegs)
    rstart = 0
    rend = rstart + figrowlength

    for n in range(nsegs):
        if rend > alignment_width:
            rend = alignment_width
        axes = f.add_subplot(gs[n])
        axes.set_xlim(rstart-0.5, rend-0.5)
        axes.set_ylim(0, 3.1)
        seq_count = len(alignment[:,0])
        width = 0.75
        ind = np.arange(rstart, rend)

        if typ == "nt":
            element_list = utilityFunctions.getNtColours()
            colours = utilityFunctions.getNtColours()
        elif typ == "aa":
            element_list = utilityFunctions.getAAColours()
            colours = utilityFunctions.getAAColours()
        height_list = {}

        for element in element_list:
            height_list[element] = []

            # A_height = []
            # G_height = []
            # C_height = []
            # U_height = []

        bottom_height = []

        for i in range(rstart, rend):
            unique, counts = np.unique(alignment[:,i], return_counts=True)
            count = dict(zip(unique, counts))
            height_per_base, info_per_base = calc_entropy(count, seq_count, that_letter, typ)
            bottom_height.append(0)

            for base, height in height_per_base.items():
                height_list[base].append(height_per_base[base])

            # A_height.append(height_per_base["A"])
            # G_height.append(height_per_base["G"])
            # C_height.append(height_per_base["C"])
            # U_height.append(height_per_base[that_letter])

            
        # plt.bar(ind, A_height, width, color=colours['A'])
        # # for idx,rect in enumerate(bar_plot):
        # #     height = rect.get_height()
        # #     ax.text(rect.get_x() + rect.get_width(), height, "A", ha='center', va='bottom')
        #
        # plt.bar(ind, G_height, width, bottom=A_height, color=colours['G'])
        # plt.bar(ind, C_height, bottom=[i+j for i,j in zip(A_height, G_height)], width=width, color=colours['C'])
        # plt.bar(ind, U_height, bottom=[i+j+k for i,j,k in zip(A_height, G_height, C_height)], width=width,
        #                                color=colours['U'])

        for base, height in height_list.items():
                plt.bar(ind, height, width, bottom=bottom_height, color=colours[base])
                bottom_height = [i+j for i,j in zip(bottom_height,height)]


        plt.xticks([rstart, rend-1], [rstart+1, rend])
        plt.yticks(np.arange(0, 3.1, 1))
        plt.xlabel("Position")
        plt.ylabel("Bit Score")

        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        rstart += figrowlength
        rend += figrowlength
    plt.savefig(figname, bbox_inches='tight', dpi=figdpi)
    plt.close()
    #todo tidy up and use normal names and normal files



def calc_entropy(count, seq_count, that_letter, typ):

    # total number of Sequences - gap number
    # adjust total height later to make up for gaps - i think that's covered (?)
    if typ == "nt":
        element_list = utilityFunctions.getNtColours()
        s = 4
        max_entropy = log(4,2)
    elif typ == "aa":
        element_list = utilityFunctions.getAAColours()
        s = 20
        max_entropy = log(20,2)

    info_per_base = {}
    freq_per_base = {}
    height_per_base = {}
    entropy_per_base = {}

    for element in element_list:
        info_per_base[element] = 0
        freq_per_base[element] = 0
        height_per_base[element] = 0
        entropy_per_base[element] = 0


    sample_size_correction = (1/log(2)) * ((s-1)/(2*seq_count)) # wahh here or later?
    gap_correction = seq_count
    if count.get("-"):
        seq_count -= count.get("-")
    gap_correction = seq_count/gap_correction
    # print('gap correction', gap_correction)
    # info_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    # freq_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    # height_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    # entropy_per_base = {"A": 0, "G": 0, that_letter: 0, "C": 0}
    entropy = 0
    if seq_count == 0:
        return height_per_base, info_per_base

    for base, quantity in count.items():
        if base != "-":
            frequency = quantity/seq_count
            freq_per_base[base] = frequency
            entropy -= frequency*log(frequency, 2)
            info_per_base[base] = max_entropy + frequency*log(frequency, 2)
            entropy_per_base[base] = -frequency*log(frequency,2)
    information_per_column = max_entropy-entropy-sample_size_correction
    # print("info", info_per_base)
    for base, quantity in info_per_base.items():
        if freq_per_base[base]*information_per_column < 0:
            height_per_base[base] = 0
        else:
            #scale to accomodate gaps, or do we even need to that here if already
            #used in the sample size correction?
            height_per_base[base] = gap_correction*freq_per_base[base]*information_per_column

    return height_per_base, info_per_base
