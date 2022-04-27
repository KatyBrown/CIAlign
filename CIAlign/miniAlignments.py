#! /usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
try:
    import CIAlign.utilityFunctions as utilityFunctions
except ImportError:
    import utilityFunctions
import math
matplotlib.use('Agg')


def arrNumeric(arr, typ):
    '''
    Converts the sequence array into a numerical matrix and a colour map
    which matplotlib can interpret as an image (similar to
                                                https://bit.ly/2CIKOEr)
    The rows in the array are inverted so that the output image has the rows
    in the same order as the input alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    Returns
    -------
    arr2: np.array
        The flipped alignment as an array of integers
    cmap: matplotlib.colors.ListedColormap
        A colour map with the colours corresponding to each base
        or amino acid
    '''
    # turn the array upside down
    arr = np.flip(arr, axis=0)
    if typ == 'nt':
        D = utilityFunctions.getNtColours()
    else:
        D = utilityFunctions.getAAColours()

    # retrieve the colours for the colour map
    keys = list(D.keys())
    ali_height, ali_width = np.shape(arr)

    # make a dictionary where each integer corresponds to a base or nt
    i = 0
    nD = dict()
    colours = []
    for key in keys:
        if key in arr:
            nD[key] = i
            colours.append(D[key])
            i += 1

    arr2 = np.empty([ali_height, ali_width])

    for x in range(ali_width):
        for y in range(ali_height):
            # numeric version of the alignment array
            arr2[y, x] = nD[arr[y, x]]

    cmap = matplotlib.colors.ListedColormap(colours)
    return (arr2, cmap)


def drawMarkUp(a, markupdict, nams, ali_width, ali_height):
    '''
    Add the "markup" to the mini alignment - on the input alignment image
    use coloured lines to show which rows, columns and positions have
    been removed by each of the parsing functions.

    Colours were selected because they are not too similar to any of
    the nucleotide or amino acid colours.

    Parameters
    ----------
    a: matplotlib.pyplot.subplot
        An open subplot containing the mini alignment for the input file
    markupdict: dict
        A dictionary where the keys are function names and the values are: for
        functions which remove whole columns a list of absolute values of
        column numbers which were removed, for functions which remove whole
        rows a list of sequence names which were removed, for functions which
        remove positions a dictionary where keys are row names and values are
        lists of absolute values of positions which were removed from that
        row
    nams: list
        The names of the sequences in the alignment
    ali_width: int
        The number of columns in the input alignment
    ali_height: int
        The number of rows (sequences) in the input alignment

    Returns
    -------

    '''
    colD = utilityFunctions.getMarkupColours()
    lineweight_h = 5 / ali_height
    # removes single positions
    if "crop_ends" in markupdict:
        colour = colD['crop_ends']
        for nam, boundary in markupdict['crop_ends'].items():
            y = len(nams) - nams.index(nam) - 1
            # left end of the sequence
            if boundary[0].shape[0] > 0:
                a.add_patch(
                        matplotlib.patches.Rectangle(
                                [boundary[0][0]-0.5, y-0.5],
                                (boundary[0][-1] - boundary[0][0]) + 1,
                                1, color=colour, lw=0, zorder=50))
                a.hlines(y, boundary[0][0]-0.5, boundary[0][-1]+0.5,
                         zorder=51, color='black',
                         lw=lineweight_h)

            # right end of the sequence
            if boundary[1].shape[0] > 0:
                a.add_patch(
                        matplotlib.patches.Rectangle(
                                [boundary[1][0]-0.5, y-0.5],
                                (boundary[1][-1] - boundary[1][0]) + 1,
                                1, color=colour, lw=0, zorder=50))
                a.hlines(y, boundary[1][0]-0.5, boundary[1][-1]+0.5,
                         zorder=51, color='black',
                         lw=lineweight_h)

    # removes whole rows
    if "remove_divergent" in markupdict:
        colour = colD['remove_divergent']
        for row in markupdict['remove_divergent']:
            y = len(nams) - nams.index(row) - 1.5
            a.add_patch(matplotlib.patches.Rectangle(
                    (-0.5, y), ali_width, 1,
                    color=colour, zorder=48, lw=0))
            a.hlines(y+0.5, -0.5, ali_width-0.5, zorder=49, color='black',
                     lw=lineweight_h)

    # removes whole columns
    if "remove_insertions" in markupdict:
        colour = colD['remove_insertions']

        for col in markupdict['remove_insertions']:

            a.add_patch(matplotlib.patches.Rectangle(
                    (col-0.5, -0.5), 1, ali_height, color=colour, zorder=46,
                    lw=0))
            for row in np.arange(ali_height):
                a.hlines(row, col-0.5, col+0.5, zorder=47, color='black',
                         lw=lineweight_h)

    # removes whole rows
    if "remove_short" in markupdict:
        colour = colD['remove_short']
        for row in markupdict['remove_short']:
            y = len(nams) - nams.index(row) - 1.5
            a.add_patch(matplotlib.patches.Rectangle(
                    (-0.5, y), ali_width, 1,
                    color=colour, zorder=44, lw=0))
            a.hlines(y+0.5, -0.5, ali_width-0.5, zorder=45, color='black',
                     lw=lineweight_h)
    # removes whole columns
    if "remove_gaponly" in markupdict:
        colour = colD['remove_gaponly']
        for col in markupdict['remove_gaponly']:
            a.add_patch(matplotlib.patches.Rectangle((col-0.5, -0.5), 1,
                                                     ali_height,
                                                     color=colour,
                                                     zorder=42, lw=0))
            for row in np.arange(ali_height):
                a.hlines(row, col-0.5, col+0.5, zorder=43, color='black',
                         lw=lineweight_h)


def drawMarkUpLegend(outfile):
    '''
    Draws a small legend (in a seperate image) showing which parsing step
    each colour represents in the marked up mini alignment.

    At the moment this is not parameterised - it will draw exactly the same
    thing every time.

    Parameters
    ----------
    outfile: str
        Path to the output file

    Returns
    -------
    None
    '''
    legend = plt.figure(figsize=(2, 2), dpi=100)
    leg = legend.add_subplot(1, 1, 1)
    colours = utilityFunctions.getMarkupColours()

    functions = {'crop_ends': 'Cropped Ends',
                 'remove_divergent': 'Too Divergent',
                 'remove_insertions': 'Insertions',
                 'remove_short': 'Too Short',
                 'remove_gaponly': 'Gap Only'}
    for i, (func, txt) in enumerate(functions.items()):
        leg.plot(1, 5-i, marker='.', color=colours[func], markersize=20)
        leg.text(2, 5-i, txt)
    leg.set_xlim(0.5, 3)
    leg.set_ylim(-1, 6)
    leg.set_axis_off()
    legend.gca().set_axis_off()
    leg.margins(0, 0)
    legend.savefig("%s_legend.png" % (outfile),
                   dpi=100, bbox_inches='tight')


def drawMiniAlignment(arr, nams, log, outfile, typ,
                      dpi=300, title=None, width=5, height=3, markup=False,
                      markupdict=None, ret=False, orig_nams=[],
                      keep_numbers=False, force_numbers=False):
    '''
    Draws a "mini alignment" image showing a small representation of the
    whole alignment so that gaps and poorly aligned regions are visible.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    log: logging.Logger
        The open log file object
    outfile: str
        Path to the output file
    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide
    dpi: int
        DPI for the output image
    title: str
        Title for the output image
    width: int
        Width of the output image
    height: int
        Height of the output image
    markup: bool
        Should the deleted rows and columns be marked on the output image?
    markupdict: dict
        Dictionary where the keys are function names and the values are
        lists of columns, rows or positions which have been removed
    ret: bool
        Return the subplot as a matplotlib object, used to make plots when
        using this function directly rather than the CIAlign workflow
    orig_nams: list
        List of names in the original input plot, used if keep_numbers is
        switched on to keep the original numbering scheme
    keep_numbers: bool
        Number the sequences (rows) based on the original CIAlign input rather
        than renumbering.

    Returns
    -------
    None

    '''
    ali_height, ali_width = np.shape(arr)

    # font size needs to scale with DPI
    fontsize = 1500 / dpi

    # what is the order of magnitude of the number of sequences (rows)
    # 0 = 1 - 10 sequences - label every row
    # 1 = 10 - 100 sequences - label every 10th row
    # 2+ = 100+ sequences - label every 100th row

    om = math.floor(math.log10(ali_height))
    tickint = 1 if om == 0 or force_numbers else 10 if om == 1 else 100

    # use thinner lines for bigger alignments
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width

    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot(1, 1, 1)
    a.set_xlim(-0.5, ali_width)
    a.set_ylim(-0.5, ali_height-0.5)

    # generate the numeric version of the array
    arr2, cm = arrNumeric(arr, typ)
    # display it on the axis
    a.imshow(arr2, cmap=cm, aspect='auto', interpolation='nearest')

    # these are white lines between the bases in the alignment - as the
    # image is actually solid
    a.hlines(np.arange(-0.5, ali_height), -0.5,
             ali_width, lw=lineweight_h, color='white',
             zorder=100)
    a.vlines(np.arange(-0.5, ali_width), -0.5,
             ali_height, lw=lineweight_v, color='white',
             zorder=100)

    # aesthetics
    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)

    if title:
        f.suptitle(title)
    for t in a.get_xticklabels():
        t.set_fontsize(fontsize)
    a.set_yticks(np.arange(ali_height-1, -1, -tickint))
    x = 1
    if tickint == 1:
        if keep_numbers:
            labs = []
            for nam in orig_nams:
                if nam in nams:
                    labs.append(x)
                x += 1
            a.set_yticklabels(labs,
                              fontsize=fontsize)
        else:
            a.set_yticklabels(np.arange(1, ali_height+1, tickint),
                              fontsize=fontsize)
    else:
        a.set_yticklabels(np.arange(0, ali_height, tickint), fontsize=fontsize)

    if markup:
        a = drawMarkUp(a, markupdict, nams, ali_width, ali_height)
        drawMarkUpLegend(outfile.replace(".png", ""))
    f.tight_layout()
    f.savefig(outfile, dpi=dpi, bbox_inches='tight')
    if ret:
        return (f)
    plt.close()
