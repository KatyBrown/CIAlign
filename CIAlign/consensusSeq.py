#! /usr/bin/env python

import matplotlib
from math import log
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import operator
import matplotlib.patheffects
import math
from matplotlib import gridspec
import copy
try:
    import CIAlign.utilityFunctions as utilityFunctions
except ImportError:
    import utilityFunctions
import os
import scipy.stats
import copy
import operator
import pandas as pd
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


def getLetters(typ='nt', fontname='monospace', dpi=500, palette="CBS"):
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
        colours = utilityFunctions.getNtColours(palette)
    elif typ == 'aa':
        colours = utilityFunctions.getAAColours(palette)
    # for each possible base/aa create temporary plot
    for base in colours.keys():
        f = plt.figure(figsize=(0.8, 1), dpi=dpi, edgecolor='black')
        a = f.add_subplot(1, 1, 1)
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        fs = getFontSize(f, a, 1)
        if base != 'Q':
            a.text(0.5, 0.02, base, fontsize=fs*0.95,
                   fontdict={'family': 'monospace',
                             'name': fontname},
                   color=colours[base], va='baseline',
                   ha='center')
        else:
            a.text(0.5, 0.1, base, fontsize=fs*0.85,
                   fontdict={'family': 'monospace',
                             'name': fontname},
                   color=colours[base], va='baseline',
                   ha='center')
        a.set_ylim(0, 1)
        plt.gca().set_axis_off()
        a.margins(0, 0)
        f.subplots_adjust(top=1, bottom=0, right=1, left=0,
                          wspace=None, hspace=0)
        a.set_frame_on(False)
        # temporarily save plot in working directory
        base = base.replace("*", "stop")
        f.savefig("%s_temp.png" % base, dpi=500,
                  pad_inches=0.5)
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


def makeLinePlot(stat, dest, ylab, dpi=300, height=3, width=5,
                 colour='#007bf5'):
    '''
    Creates a plot of the coverage

    Parameters
    ----------
    coverage: list of strings
        Coverage for alignment

    dest: str
        folder to store file

    ylab: str
        label for y axis

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
    x = np.arange(0, len(stat), 1)
    y = stat

    xmax = x.max()
    ymax = max(y)

    x_div = 10**(math.floor(np.log10(xmax))) / 5
    y_div = 10**(math.floor(np.log10(ymax)))

    # used for polynomial interpolation
    # xmin = x.min()
    # N = 100
    # xx = np.linspace(xmin, xmax, N)

    # plain plotting of the coverage
    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot(2, 1, 1)
    a.plot(x, y, color=colour, lw=1)
    a.set_xlabel('Position', fontsize=fontsize)
    a.set_ylabel(ylab, fontsize=fontsize)

    b = f.add_subplot(2, 1, 2)

    # 30 represents number of points to make between T.min and T.max
    # interpolating the coverage function to make it smooth
    xnew = np.linspace(x.min(), x.max(), 30)
    spline = interpolate.make_interp_spline(x, y, k=3)
    b.plot(xnew, spline(xnew), color=colour, lw=1)
    b.set_xlabel('Position', fontsize=fontsize)
    b.set_ylabel('%s (Smoothed)' % ylab, fontsize=fontsize)

    # polynomial interpolation leaving this in just in case
    # c = f.add_subplot('313')
    # z = np.polyfit(x, bla, 30)
    # p = np.poly1d(z)
    # c.plot(xx, p(xx))
    for sp in f.axes:
        sp.set_xticks(np.arange(0, xmax*1.1, x_div))
        sp.set_xticklabels([int(x) for x in np.arange(0, xmax*1.1, x_div)],
                           fontsize=fontsize,
                           rotation='vertical')
        sp.set_yticks(np.arange(0, ymax*1.1, y_div))
        sp.set_yticklabels([
            "%.1f" % y for y in np.arange(0, ymax*1.1, y_div)],
            fontsize=fontsize)
        sp.set_xlim(0, xmax)
    f.tight_layout()
    f.savefig(dest, dpi=dpi, bbox_inches='tight')


def sequence_logo(alignment,
                  figname,
                  typ='nt',
                  figfontname='Arial',
                  figdpi=300,
                  figrowlength=50,
                  start=0,
                  end=0,
                  palette='CBS',
                  return_fig=False):
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

    start: int
           start pos to be turned into logo

    end: int
         end pos to be turned into logo

    Returns
    -------
    none
    '''

    if start == 0 and end == 0:
        alignment_width = len(alignment[0, :])
    else:
        if end == 0:
            end = len(alignment[0, :])
        alignment_width = len(alignment[0, start:end])

    if alignment_width < figrowlength:
        figrowlength = alignment_width
    nsegs = math.ceil(alignment_width / figrowlength)
    f = plt.figure(figsize=(figrowlength, nsegs*4), dpi=figdpi)
    gs = gridspec.GridSpec(ncols=1, nrows=nsegs)
    getLetters(typ=typ, fontname=figfontname, dpi=figdpi, palette=palette)
    rstart = start
    rend = rstart + figrowlength
    for n in range(nsegs):

        if rend > (alignment_width + start):
            rend = alignment_width + start
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
            height_per_base, info_per_base, fq = calc_entropy(
                count, len(alignment[:, 0]), typ=typ)
            height_sum_higher = 0
            Z = zip(height_per_base.keys(), height_per_base.values())
            Z = sorted(Z, key=lambda x: x[1])
            for base, height in Z:
                if height > 0:
                    b = base.replace("*", "stop")
                    L = plt.imread("%s_temp.png" % b)
                    a.imshow(L, extent=(i, i+1, height_sum_higher,
                                        height_sum_higher+height),
                             filternorm=False)

                    height_sum_higher += height
        a.axis(limits)
        if rend - start < 25:
            a.set_xticks([int(x) for x in np.arange(rstart-0.5, rend, 1)])
            a.set_xticklabels([int(x) for x in np.arange(rstart, rend+0.5, 1)])
        else:
            a.set_xticks([rstart, rend])
            a.set_xticklabels([rstart, rend])
        a.set_xlim(rstart, rstart+figrowlength)
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        if n == (nsegs - 1):
            a.set_xlabel("Position")
        a.set_ylabel("Bit Score")
        rstart += figrowlength
        rend += figrowlength

    # obtain colours
    if typ == 'nt':
        allbases = utilityFunctions.getNtColours(palette=palette)
    elif typ == 'aa':
        allbases = utilityFunctions.getAAColours(palette=palette)
    for base in allbases:
        b = base.replace("*", "stop")
        os.unlink("%s_temp.png" % b)
    # save plot using figname
    if return_fig:
        return (f)
    f.savefig(figname, dpi=figdpi, bbox_inches='tight')
    plt.close()


def sequence_bar_logo(alignment,
                      figname,
                      typ='nt',
                      figdpi=300,
                      figrowlength=50,
                      start=0,
                      end=0,
                      palette='CBS',
                      return_fig=False):
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

    start: int
           start pos to be turned into logo

    end: int
         end pos to be turned into logo

    Returns
    -------
    none
    '''

    if start == 0 and end == 0:
        alignment_width = len(alignment[0, :])
    else:
        if end == 0:
            end = len(alignment[0, :])
        alignment_width = len(alignment[0, start:end])

    if alignment_width < figrowlength:
        figrowlength = alignment_width
    nsegs = math.ceil(alignment_width / figrowlength)
    f = plt.figure(figsize=(figrowlength/5, nsegs*2), dpi=figdpi)
    gs = gridspec.GridSpec(ncols=1, nrows=nsegs)
    rstart = start
    rend = rstart + figrowlength

    for n in range(nsegs):
        if rend > (alignment_width + start):
            rend = alignment_width + start
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
            element_list = utilityFunctions.getNtColours(palette=palette)
            colours = utilityFunctions.getNtColours(palette=palette)
        elif typ == "aa":
            element_list = utilityFunctions.getAAColours(palette=palette)
            colours = utilityFunctions.getAAColours(palette=palette)
        height_list = {}
        for element in element_list:
            height_list[element] = []

        bottom_height = []
        # for each column calculate heights via entropy
        # and scale letters accordlingly
        for i in range(rstart, rend):
            unique, counts = np.unique(alignment[:, i], return_counts=True)
            count = dict(zip(unique, counts))
            height_per_base, info_per_base, fr = calc_entropy(count, seq_count,
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
    if return_fig:
        return (f)
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
    if seq_count == 0:
        return height_per_base, info_per_base, 0
    # correct for small sample sizes

    sample_size_correction = (s-1) / (2 * np.log(2) * seq_count)
    gap_correction = seq_count
    if count.get("-"):
        seq_count -= count.get("-")

    # correct for gaps, since they lower the information content
    gap_correction = seq_count/gap_correction

    entropy = 0

    # calculate entropy, from that information, from that height
    freqs = []
    for base, quantity in count.items():
        if base != "-":
            frequency = quantity/seq_count
            freq_per_base[base] = frequency
            entropy -= frequency*log(frequency, 2)
            info_per_base[base] = max_entropy + frequency*log(frequency, 2)
            entropy_per_base[base] = -frequency * log(frequency, 2)
            freqs.append(frequency)

    information_per_column = max_entropy - entropy - sample_size_correction

    # if the information content is constant throughout the column,
    # these value will be negative. Since this does not add any information
    # set them to 0
    # they can be negative due to the sample size correction
    # (otherwise they'd be 0)
    for base, quantity in info_per_base.items():
        if freq_per_base[base]*information_per_column < 0:
            height_per_base[base] = 0
        else:
            # scale to accomodate gaps
            height_per_base[base] = (gap_correction *
                                     freq_per_base[base] *
                                     information_per_column)

    return height_per_base, info_per_base, freqs


def calcConservationAli(alignment, typ):
    '''
    Calculate alignment conservation as the heights the letters would
    be in a sequence logo.

    Parameters
    -----------
    alignment: np.array
        The alignment stored as a numpy array

    typ: str
        nt or aa

    Returns
    -------
    heights: list
        A list containing the letter height for each column
    ents: list
        A list containing the entropy calculated for each column
    '''
    alignment_width = len(alignment[0, :])
    seq_count = len(alignment[:, 0])
    heights_per_col = []
    ents = []
    for i in range(0, alignment_width):
        unique, counts = np.unique(alignment[:, i], return_counts=True)
        count = dict(zip(unique, counts))
        height_per_base, info_per_base, freqs = calc_entropy(count, seq_count,
                                                             typ)
        heights_per_col.append(height_per_base)
        ent = scipy.stats.entropy(freqs)
        ents.append(ent)
    heights = [sum(x.values()) for x in heights_per_col]
    return (heights, ents)


def compareAlignmentConsensus(arr, typ, booleanOrSimilarity="boolean",
                              MatrixName="default"):
    '''
    Compares the alignment of the input alignment to the consensus of that
    array, and will either output a boolean array, or will use a similarity
    matrix, specified as MatrixName to output an array containing the
    similarity scores of the alignment compared to the consensus.

    Parameters
    ----------
    arr: np.array
        the aligned sequences

    typ: str
        nt or aa

    booleanOrSimilarity: str
        Create a boolean or similarity matrix (default = 'boolean')

    MatrixName: str
        substitution matrix name (default = 'default')

    Returns
    -------
    new_arr: np.array
        A boolean array showing if positions in the alignment are identical
        to the consensus

    new_Sarr: np.array
        A integer array containing the similarity scores for the
        alignment compared to the consensus calculated via a scoring matrix.
    '''
    consensus, _ = np.array(findConsensus(arr, '',
                                          consensus_type='majority_nongap'))
    ci_dir = os.path.dirname(utilityFunctions.__file__)
    matrix_dir = "%s/similarity_matrices" % (ci_dir)

    if booleanOrSimilarity == "boolean":
        bool_array = np.array([])
        bool_arrL = np.empty(dtype=bool, shape=(0, len(consensus)))
        # declares the numpy arrays
        for e in range(1, (len(arr[:, 0])+1)):
            # iterates over the rows of the sequences
            z = e - 1
            for i in range(1, (len(arr[0, :])+1)):
                # iterates over the columns of the sequences
                x = i-1
                if arr[z, x] == consensus[x]:
                    # verifies if the current value being iterated is equal to
                    # the equivalent value inline with the consensus
                    bool_array = np.append(bool_array, [True], axis=None)
                elif arr[z, x] == "-" or consensus[x] == "-":
                    bool_array = np.append(bool_array, [float('nan')],
                                           axis=None)
                else:
                    bool_array = np.append(bool_array, [False], axis=None)
            bool_arrL = np.vstack([bool_arrL, bool_array])
            bool_array = np.array([])

        new_arr = copy.deepcopy(bool_arrL)
        new_arr = bool_arrL.astype(float)
        # returns the new boolean array containing the verified alignment to
        # the consensus
        return new_arr
    else:
        # generates the consensus
        Sarray = np.array([])
        SarrL = np.empty(dtype=int, shape=(0, len(consensus)))
        # declares the numpy arrays
        tab = pd.read_csv("%s/matrices.txt" % ci_dir, sep="\t", index_col=0)
        if typ == "aa":
            # verifies if the typ is amino acid or nucleotide
            if MatrixName != "default":
                if tab.loc[MatrixName][0] != typ:
                    raise RuntimeError("This substitution matrix is not \
                                       valid for \
                                       alignment type %s" % typ)
                else:
                    # verifies if the user would like to use the default matrix or
                    # their own
                    mat = pd.read_csv("%s/%s" % (matrix_dir, MatrixName),
                                      comment="#", sep=r"\s+")
            elif MatrixName == "default":
                mat = pd.read_csv("%s/BLOSUM62" % (matrix_dir),
                                  comment="#", sep=r"\s+")
        elif typ == "nt":
            if MatrixName != "default":
                if tab.loc[MatrixName][0] != typ:
                    raise RuntimeError("This substitution matrix is not \
                                       valid for alignment type %s" % typ)
                else:
                    # verifies if the user would like to use the default
                    # matrix or their own
                    mat = pd.read_csv("%s/%s" % (matrix_dir, MatrixName),
                                      comment="#", sep=r"\s+")
            elif MatrixName == "default":
                mat = pd.read_csv("%s/NUC.4.4" % (matrix_dir), comment="#",
                                  sep=r"\s+")
        for e in range(1, (len(arr[:, 0])+1)):
            # iterates over the rows of the sequences
            z = e-1
            for i in range(1, (len(arr[0, :])+1)):
                #  iterates over the columns of the sequences
                x = i-1
                if not arr[z, x] == "-":
                    if arr[z, x] == "U" and typ == 'nt':
                        thischar = "T"
                    else:
                        thischar = arr[z, x]
                    if consensus[x] == "U" and typ == 'nt':
                        conschar = "T"
                    else:
                        conschar = consensus[x]
                    score = mat.loc[thischar, conschar]
                    Sarray = np.append(Sarray, [score])
                elif arr[z, x] == "-":
                    # sets the value of '-' as 0
                    Sarray = np.append(Sarray, float('nan'))
            SarrL = np.vstack([SarrL, Sarray])
            Sarray = np.array([])
        new_Sarr = copy.deepcopy(SarrL)
        new_Sarr = SarrL.astype(float)
        # returns the new similarity array containing the verified alignment
        # to the consensus
        return new_Sarr


def plotResidueFrequencies(arr, typ, outfile, dpi=300, width=3, height=5):
    '''
    Plots a bar chart of the percentage of each base/aa in the
    alignment, also plots a graph comparing the
    proportion of C/G bases to A/T bases.

    Parameters
    -----------
    arr: np.array
        The sequence stored as a numpy array

    typ: str
        nt or aa

    outfile: str
        Path to the file to store the graph image

    dpi: int
        Dots per inch for output file

    width: int
        Output image width in inches. For amino acids this will be x3
        (to allow same default for both nt and aa).

    height: int
        Output image height in inches.

    Outputs
    -------
    The bar chart in the file "outfile""

    '''
    # Get the residues and colours
    if typ == "nt":
        Slist = list(utilityFunctions.getNtColours().keys())[0:5]
        Scols = list(utilityFunctions.getNtColours().values())[0:5]
    elif typ == "aa":
        Slist = list(utilityFunctions.getAAColours().keys())[0:20]
        Scols = list(utilityFunctions.getAAColours().values())[0:20]

    # Set up empty arrays
    baseT = np.array([])
    CGcount = float(0)
    ATcount = float(0)

    # This allows the plots to be on a similar scale with the same default
    # width even though there are more amino acids
    if typ == 'aa':
        width = width * 3

    # Create the blank figure, add space between plots
    f = plt.figure(figsize=(width, height))
    plt.subplots_adjust(hspace=0.5)

    for i in range(0, len(Slist)):
        # Count C and G vs A and T
        if Slist[i] == 'C' or Slist[i] == 'G':
            CGcount = CGcount + np.sum(arr == Slist[i])
        if Slist[i] == 'A' or Slist[i] == 'T' or Slist[i] == 'U':
            ATcount = ATcount + np.sum(arr == Slist[i])

        # Sub "U" for "T" for RNA
        if typ == 'nt' and Slist[i] == 'T':
            thischar = "U"
        else:
            thischar = Slist[i]

        # Calculate the proportion of this character overall
        prop = (np.sum(arr == thischar) / np.size(arr))
        baseT = np.append(baseT, prop)

    # For nt, plot counts and CG frequency
    # For aa, just plot counts
    if typ == 'nt':
        subplot1 = f.add_subplot(2, 1, 1)
        subplot1.set_xlabel("Nucleotide")
    else:
        subplot1 = f.add_subplot(1, 1, 1)
        subplot1.set_xlabel("Amino Acid")

    # Plot the bar chart
    subplot1.bar(Slist, baseT, color=Scols, width=0.4)

    # Aesthetics for plot
    subplot1.set_ylabel("Proportion")
    subplot1.spines['right'].set_visible(False)
    subplot1.spines['top'].set_visible(False)
    subplot1.set_title("Proportion of Residues in Alignment")

    # Add the AT vs CG plot for amino acids
    if typ == 'nt':
        subplot2 = f.add_subplot(2, 1, 2)
        # Differentiate U and T
        if "U" in arr:
            label = "A or U"
            subplot2.set_title("CG and AU Proportion")
        else:
            label = "A or T"
            subplot2.set_title("CG and AT Proportion")

        # Plot the bar chart
        subplot2.bar(["C or G", label], [(CGcount / (ATcount + CGcount)),
                                         (ATcount / (ATcount + CGcount))],
                     color="Orange", width=0.4)

        # Aesthetics for plot
        subplot2.spines['right'].set_visible(False)
        subplot2.spines['top'].set_visible(False)
        subplot2.set_ylabel("Proportion")
        subplot2.set_xlabel("Nucleotides")

    # Save as outfile
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()


def residueChangeCount(arr, typ, outfile, dpi=300, width=3, height=5):
    '''
    Plots a bar chart of the percentage of each possible change
    between the consensus sequence and sequences in the alignment.

    Parameters
    -----------
    arr: np.array
        The alignment stored as a numpy array

    typ: str
        Type - aa or nt. This graph is plotted for nt only

    outfile: str
        Path to the file to store the graph image

    dpi: int
        Dots per inch for output file

    width: int
        Output image width in inches

    height: int
        Output image height in inches

    Outputs
    -------
    The bar chart in the file "outfile""
    '''

    # Get the consensus
    consensus, _ = np.array(findConsensus(
        arr, '', consensus_type='majority_nongap'))

    # Track all the differences between the consensus and the array
    # (excluding gaps)
    NDiff = np.array([])
    for i in range(1, len(arr[:, 0])):
        z = i - 1
        for e in range(1, len(arr[0, :])):
            x = e - 1
            if arr[z, x] != consensus[x] and consensus[x] != "-":
                if arr[z, x] != "-":
                    NDiff = np.append(NDiff, str(arr[z, x] + consensus[x]))
    total_arr = np.size(NDiff)
    xAx = ["AT", "AC", "AG", "TA", "TC", "TG",
           "CA", "CT", "CG", "GA", "GC", "GT"]
    yAy = np.array([])
    for o in range(1, len(xAx)+1):
        w = o-1
        yAy = np.append(yAy, (sum((NDiff == str(xAx[w]))[:])/total_arr))

    f = plt.figure(figsize=(width*3, height))
    a = f.add_subplot(1, 1, 1)
    a.bar(xAx, yAy, color="orange", width=0.4)

    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.set_ylabel("Proportion")
    a.set_xlabel("Change (Consensus to Alignment)")

    # Save as outfile
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()
