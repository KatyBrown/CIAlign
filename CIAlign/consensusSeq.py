#! /usr/bin/env python

import matplotlib
from math import log
import numpy as np
# import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import operator
import matplotlib.patheffects
import math
from matplotlib import gridspec
try:
    import CIAlign.utilityFunctions as utilityFunctions
except ImportError:
    import utilityFunctions
import os
matplotlib.use('Agg')


def getAxisUnits(subplot):
    '''
    Translates the height of one unit of y axis
    and width of one unit of x axis into pixels
    Translates height of y axis and width of x axis
    to pixels

    Parameters
    ----------
    subplot: matplotlib.pyplot.subplot
        An open subplot

    Returns
    -------
    D: dict
        Keys are names of axis heights and widths, the values are in units
        or pixels
    '''
    axis_dimensions = subplot.transData.transform(
        [(subplot.get_xlim()[1],
          subplot.get_ylim()[1]), (0, 0)]) - subplot.transData.transform((
              0, 0))
    D = dict()
    # height of y in pixels
    D['axis_height_px'] = axis_dimensions[0][1]
    # height of x in pixels
    D['axis_width_px'] = axis_dimensions[0][0]
    D['axis_bottom'], D['axis_top'] = subplot.get_ylim()
    D['axis_left'], D['axis_right'] = subplot.get_xlim()
    # total number of units on x axis
    D['axis_width_u'] = D['axis_right'] - D['axis_left']
    # total number of units on y axis
    D['axis_height_u'] = D['axis_top'] - D['axis_bottom']
    # height of one y axis unit in pixels
    D['u_height_px'] = D['axis_height_px'] / D['axis_height_u']
    # width of one x axis unit in pixels
    D['u_width_px'] = D['axis_width_px'] / D['axis_width_u']

    return (D)


def getFontSize(figure, subplot, height_u):
    '''
    Based on axes units, converts specified number of y axis units into
    the font size in points

    Parameters
    ----------
    figure: matplotlib.Figure
        An open figure

    subplot: matplotlib.pyplot.subplot
        An open subplot

    height_u: float
        Height in units

    Returns
    -------
    Height_points*1.8: float
        adjusted height in points
    '''

    D = getAxisUnits(subplot)
    perc_height = height_u / (D['axis_top'] - D['axis_bottom'])
    height_pixels = D['axis_height_px'] * perc_height
    dpi = figure.dpi
    height_points = PixelsToPoints(height_pixels, dpi)

    return (height_points * 1.8)


def PixelsToPoints(pixels, dpi):
    '''
    Converts pixel values to point values

    Parameters
    ----------
    pixels: float
        pixel value

    dpi: int
        dpi value

    Returns
    -------
    points: float
        points value
    '''
    # one point is 1/72 of an inch
    points = pixels * (72 / dpi)

    return (points)


def getLetters(typ='nt', fontname='monospace', dpi=500):
    '''
    Generates a temporary image file for every letter (4 nt or 20 aa)
    Each letter extends to the full length of both axes

    Parameters
    ----------
    typ: string
        nt (default) or aa

    fontname: string
        name of a font (default: monospace)

    dpi: int
        DPI value (default: 500)

    Returns
    -------
    none
    '''
    # obtain color scheme depending on nt or aa alignment
    if typ == 'nt':
        colours = utilityFunctions.getNtColours()
    elif typ == 'aa':
        colours = utilityFunctions.getAAColours()
    # for each possible base/aa create temporary plot
    for base in colours.keys():
        f = plt.figure(figsize=(1, 1), dpi=dpi, edgecolor='black')
        a = f.add_subplot(111)
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        fs = getFontSize(f, a, 1)
        a.text(0.5, 0.02, base, fontsize=fs*0.95, fontdict={'family': 'monospace',
                                                    'name': fontname},
               color=colours[base], va='baseline', ha='center')
        plt.gca().set_axis_off()
        a.margins(0, 0)
        f.subplots_adjust(top=1, bottom=0, right=1, left=0,
                          wspace=None, hspace=None)
        a.set_frame_on(False)
        # temporarily safe plot in working directory
        f.savefig("%s_temp.png" % base, dpi=500,
                  pad_inches=0.1)
        plt.close()


def findConsensus(alignment, log, consensus_type="majority"):
    '''
    Calculates the consensus sequence and the coverage for the alignment
    Utilises different types of the consensus

    Parameters
    ----------
    alignment: np.array
        The alignment stored as a numpy array

    log: string
        name of log file

    consensus_type: string
        majority (default) or majority_nongap

    Returns
    -------
    consensus: list of strings
        consensus sequence

    coverage: list of floats
        alignment coverage
    '''

    consensus = []
    coverage = []
    numberOfSequences = len(alignment[:, 0])

    # need the reverse of array to access every column
    for i in range(0, len(alignment[0, :])):
        unique, counts = np.unique(alignment[:, i], return_counts=True)
        count = dict(zip(unique, counts))
        unique_ng = unique[unique != "-"]
        counts_ng = counts[unique != "-"]
        # deal with gap only columns
        if counts_ng.size == 0:
            count_ng = {"N": len(alignment[:, i])}
            nonGapContent = 0
        else:
            count_ng = dict(zip(unique_ng, counts_ng))
            if '-' in count:
                nonGapContent = 1-(count['-']/numberOfSequences)
            else:
                nonGapContent = 1

        # dealing with gap only collumns
        maxChar, maxCount = max(count.items(), key=operator.itemgetter(1))
        maxChar_ng, maxCount_ng = max(count_ng.items(),
                                      key=operator.itemgetter(1))

        # if there is an equal number of gap and non-gap characters at the
        # site, keep the non-gap character
        # if majoriy_nongap chosen, use the nongap
        if maxCount_ng == maxCount or consensus_type == "majority_nongap":
            maxChar = maxChar_ng
        consensus.append(maxChar)
        coverage.append(nonGapContent)

    return consensus, coverage


def makeCoveragePlot(coverage, dest, dpi=300, height=3, width=5,
                     colour='#007bf5'):
    '''
    Creates a plot of the coverage

    Parameters
    ----------
    coverage: list of strings
        Coverage for alignment

    dest: str
        folder to store file

    dpi: int
        DPI value (default: 500)

    height: int
            height of plot, default: 3

    width: int
            height of plot, default: 5

    colour: str
            coverage colour (default: #007bf5)

    Returns
    -------
    none
    '''

    fontsize = 1500 / dpi
    x = np.arange(0, len(coverage), 1)
    y = coverage

    xmax = x.max()
    # used for polynomial interpolation
    # xmin = x.min()
    # N = 100
    # xx = np.linspace(xmin, xmax, N)

    # plain plotting of the coverage
    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot('211')
    a.plot(x, y, color=colour)
    a.set_xlabel('Position', fontsize=fontsize)
    a.set_ylabel('Coverage', fontsize=fontsize)
    a.set_xticks([0, xmax])
    a.set_xticklabels([0, xmax], fontsize=fontsize)
    a.set_yticks(np.arange(0, 1.1, 0.5))
    a.set_yticklabels(np.arange(0, 1.1, 0.5), fontsize=fontsize)

    # b = f.add_subplot('212')

    # polynomial interpolation leaving this in just in case
    # c = f.add_subplot('313')
    # z = np.polyfit(x, bla, 30)
    # p = np.poly1d(z)
    # c.plot(xx, p(xx))
    # xnew = np.linspace(x.min(),x.max(),300)
    # 300 represents number of points to make between T.min and T.max

    # interpolating the coverage function to make it smooth
    # t, c, k = interpolate.splrep(x, y, s=0, k=4)
    # spline = interpolate.BSpline(t, c, k, extrapolate=False)
    # b.plot(xx, spline(xx), color=colour)
    # b.set_xlabel('Position', fontsize=fontsize)
    # b.set_ylabel('Coverage (Smoothed)', fontsize=fontsize)
    # b.set_xticks([0, xmax])
    # b.set_xticklabels([0, xmax], fontsize=fontsize)
    # b.set_yticks(np.arange(0, 1.1, 0.5))
    # b.set_yticklabels(np.arange(0, 1.1, 0.5), fontsize=fontsize)
    f.savefig(dest, dpi=dpi, bbox_inches='tight')


def sequence_logo(alignment,
                  figname,
                  typ='nt',
                  figfontname='Arial',
                  figdpi=300,
                  figrowlength=50):
    '''
    Creates a sequence logo based on an entropy calculation using letters
    Scales the letters according to the information content of the alignment
    Representation of the consensus sequence of the alignment

    Parameters
    ----------
    alignment: np.array
        The alignment stored as a numpy array

    figname: str
        name of figure

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    figfontname: str
            Name of font, default: Arial

    figdpi: int
            DPI (default: 300)

    figrowlength: int
            clength of figure (default: 50)

    Returns
    -------
    none
    '''
    alignment_width = len(alignment[0, :])
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
        if typ == 'nt':
            a.set_ylim(0, 2.1)
            a.set_yticks(np.arange(0, 2.1, 1))
        elif typ == 'aa':
            a.set_ylim(0, 4.6)
            a.set_yticks(np.arange(0, 4.6, 1))
        limits = a.axis()

        # for each column calculate heights via entropy
        # and scale letters accordlingly
        for i in range(rstart, rend):
            unique, counts = np.unique(alignment[:, i],
                                       return_counts=True)
            count = dict(zip(unique, counts))
            height_per_base, info_per_base = calc_entropy(count,
                                                          len(alignment[:, 0]),
                                                          typ=typ)
            height_sum_higher = 0
            Z = zip(height_per_base.keys(), height_per_base.values())
            Z = sorted(Z, key=lambda x: x[1])
            for base, height in Z:
                if height > 0:
                    L = plt.imread("%s_temp.png" % base)
                    a.imshow(L, extent=(i, i+1, height_sum_higher,
                                        height_sum_higher+height),
                             filternorm=False)

                    height_sum_higher += height
        a.axis(limits)
        a.set_xticks([rstart, rend])
        a.set_xticklabels([rstart, rend])
        
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        if n == (nsegs - 1):
            a.set_xlabel("Position")
        a.set_ylabel("Bit Score")
        rstart += figrowlength
        rend += figrowlength
    # obtain colours
    if typ == 'nt':
        allbases = utilityFunctions.getNtColours()
    elif typ == 'aa':
        allbases = utilityFunctions.getAAColours()
    for base in allbases:
        os.unlink("%s_temp.png" % base)
    # save plot using figname
    f.savefig(figname, dpi=figdpi, bbox_inches='tight')
    plt.close()


def sequence_bar_logo(alignment,
                      figname,
                      typ='nt',
                      figdpi=300,
                      figrowlength=50):
    '''
    Creates a sequence logo based on an entropy calculation using bars
    Scales the bars according to the information content of the alignment
    Representation of the consensus sequence of the alignment

    Parameters
    ----------
    alignment: np.array
        The alignment stored as a numpy array

    figname: str
        name of figure

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    figfontname: str
            Name of font, default: Arial

    figdpi: int
            DPI (default: 300)

    figrowlength: int
            clength of figure (default: 50)

    Returns
    -------
    none
    '''
    alignment_width = len(alignment[0, :])
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
        if typ == 'nt':
            axes.set_ylim(0, 2.1)
            axes.set_yticks(np.arange(0, 2.1, 1))
        elif typ == 'aa':
            axes.set_ylim(0, 4.6)
            axes.set_yticks(np.arange(0, 4.6, 1))
        seq_count = len(alignment[:, 0])
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

        bottom_height = []
        # for each column calculate heights via entropy
        # and scale letters accordlingly
        for i in range(rstart, rend):
            unique, counts = np.unique(alignment[:, i], return_counts=True)
            count = dict(zip(unique, counts))
            height_per_base, info_per_base = calc_entropy(count, seq_count,
                                                          typ)
            bottom_height.append(0)

            # need a list of each nt/aa separately to plot them as bars
            for base, height in height_per_base.items():
                height_list[base].append(height_per_base[base])

        # stag bars on top of each other
        for base, height in height_list.items():
            plt.bar(ind, height, width, bottom=bottom_height,
                    color=colours[base])
            bottom_height = [i+j for i, j in zip(bottom_height, height)]

        plt.xticks([rstart, rend-1], [rstart+1, rend])
        plt.yticks(np.arange(0, 2.1, 1))
        plt.xlabel("Position")
        plt.ylabel("Bit Score")

        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        rstart += figrowlength
        rend += figrowlength
    # save plot as figname
    plt.savefig(figname, bbox_inches='tight', dpi=figdpi)
    plt.close()


def calc_entropy(count, seq_count, typ):
    '''
    Creates a sequence logo based on an entropy calculation using bars
    Scales the bars according to the information content of the alignment
    Representation of the consensus sequence of the alignment

    Parameters
    ----------
    count: dict
        of nt/aa with counts

    seq_count: int
        number of sequences in alignment

    typ: str
        nt or aa

    Returns
    -------
    height_per_base: dictionary
        height for each nt/aa

    info_per_base: dictionary
        information content for each nt/aa

    '''
    # obtain nt/aa lists, use colour scheme for that w/o using colours here
    # just because another list of nt/aa would be obsolete
    if typ == "nt":
        element_list = utilityFunctions.getNtColours()
        s = 4
        max_entropy = log(4, 2)
    elif typ == "aa":
        element_list = utilityFunctions.getAAColours()
        s = 20
        max_entropy = log(20, 2)

    info_per_base = {}
    freq_per_base = {}
    height_per_base = {}
    entropy_per_base = {}

    for element in element_list:
        info_per_base[element] = 0
        freq_per_base[element] = 0
        height_per_base[element] = 0
        entropy_per_base[element] = 0

    # correct for small sample sizes
    sample_size_correction = (1/log(2)) * ((s-1)/(2*seq_count))
    gap_correction = seq_count
    if count.get("-"):
        seq_count -= count.get("-")
    # correct for gaps, since they lower the information content
    gap_correction = seq_count/gap_correction

    entropy = 0
    if seq_count == 0:
        return height_per_base, info_per_base

    # caluclate entropy, from that information, from that height
    for base, quantity in count.items():
        if base != "-":
            frequency = quantity/seq_count
            freq_per_base[base] = frequency
            entropy -= frequency*log(frequency, 2)
            info_per_base[base] = max_entropy + frequency*log(frequency, 2)
            entropy_per_base[base] = -frequency * log(frequency, 2)
    information_per_column = max_entropy-entropy-sample_size_correction

    # if the information content is constant throughout the column,
    # these value will be negative. Since this does not add any information
    # set them to 0
    for base, quantity in info_per_base.items():
        if freq_per_base[base]*information_per_column < 0:
            height_per_base[base] = 0
        else:
            # scale to accomodate gaps
            height_per_base[base] = (gap_correction *
                                     freq_per_base[base] *
                                     information_per_column)

    return height_per_base, info_per_base
