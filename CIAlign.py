#! /usr/bin/env python

import logging
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
import cropseq



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


def writeOutfile(outfile, arr, nams):
    out = open(outfile, "w")
    i = 0
    for i, nam in enumerate(nams):
        out.write(">%s\n%s\n" % (nam, "".join(list(arr[i]))))
    out.close()


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
            "-": '#FFFFFF'}


def seqType(arr):
    '''
    Detects if an alignment is of nucleotides or amino acids
    '''
    seq1 = arr[0]
    nucs = set(list(getNtColours().keys()))
    aas = set(list(getAAColours().keys()))
    n = 0
    a = 0
    x = 0
    for s in seq1:
        s = s.upper()
        if s in nucs:
            n += 1
        if s in aas:
            a += 1
        if s not in aas and s not in nucs:
            x += 1
    counts = n, a, x
    if n == max(counts):
        return "nt"
    elif a == max(counts):
        return "aa"
    else:
        raise RuntimeError ("Majority of positions are not known nucleotides or amino acids")


def removeInsertions(arr, log, min_size, max_size, min_flank):
    '''
    Removes insertions of size between min_size and
    max_size which have lower coverage than their flanking regions, if at least
    twice as many other sequences have the flanking regions
    but not the insertions.
    '''
    log.info("Removing insertions\n")
    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column
    sums = sum(boolarr)

    # run a sliding window along the alignment and check for regions
    # which have higher coverage at the ends of the window than in the
    # middle - store these in put_indels
    put_indels = set()
    for size in range(min_size+2, max_size+1, 2):
        for i in range(0, len(sums)+1 - size, 1):
            these_sums = sums[i:i+size]
            # take the number of non gap positions
            # for each column in this window
            ns = np.array(range(i, i+size))
            left = these_sums[0]
            right = these_sums[-1]
            # record sites with lower coverage than their flanking seqs
            # and less than mincov total coverage
            x = set(ns[(these_sums < left) &
                       (these_sums < right)])
            put_indels = put_indels | x
    put_indels = np.array(sorted(list(put_indels)))
    # for the putative indels, check if any sequence without the indel
    # flanks the indel site - if it does it confirms it is an indel
    rmpos = set()
    print (len(put_indels))
    for p in put_indels:
        has_flanks = 0
        for a in arr:
            nongaps = a == "-"
            thispos = nongaps[p]
            if thispos:
                leftsum = sum(nongaps[:p])
                rightsum = sum(nongaps[p:])
                if leftsum > min_flank and rightsum > min_flank:
                    has_flanks += 1
            if has_flanks == sums[p]:
                rmpos.add(p)
                break
    # make a list of positions to keep
    rmpos = np.array(list(rmpos))
    keeppos = np.arange(0, len(sums))
    keeppos = np.invert(np.in1d(keeppos, rmpos))
    log.info("Removing sites %s" % (", ".join([str(x) for x in rmpos])))
    arr = arr[:, keeppos]
    return (arr, set(rmpos))


def removeTooShort(arr, log, min_length, nams):
    '''
    Removes sequences with fewer than min_length non-gap positions from
    the alignment.
    '''
    if len(arr) != 0:
        arrT = arr.transpose()
        sums = sum(arrT != "-")
        arr = arr[sums > min_length]
        rmnames = np.array(nams)[sums <= min_length]
        log.info("Removing sequences %s" % (", ".join(list(rmnames))))
    else:
        rmnames = set()
    return (arr, rmnames)


def removeGapOnly(arr, log):
    if len(arr) != 0:
        sums = sum(arr == "-")
        rmpos = set(np.where(sums == len(arr[:,0]))[0])
        arr = arr[:, sums != len(arr[:,0])]
        log.info("Removing gap only sites %s" % (", ".join([str(x) for x in rmpos])))
    else:
        rmpos = set()
    return (arr, rmpos)


def cropEnds(arr, log, nams, mingap):
    newarr = []
    r = dict()
    for i, row in enumerate(arr):
        start, end = cropseq.determineStartEnd(row, mingap)
        start = max(start - 1, 0)
        end = end + 1
        newseq = "-" * start + "".join(row[start:end]) + "-" * (len(row) - end)
        newseq = np.array(list(newseq))
        s = sum(newseq != row)
        if s != 0:
            nam = nams[i]
            non_gap_start = sum(row[0:start] != "-")
            non_gap_end = sum(row[end:] != "-")
            if non_gap_start != 0:
                log.info("Removed %i bases from start of %s" % (non_gap_start, nam))
            if non_gap_end != 0:
                log.info("Removed %i bases from end of %s" % (non_gap_end, nam))
            startpos = np.where(row[0:start] != "-")[0]
            endpos = np.where(row[end:] != "-")[0] + end
            r[nam] = ((startpos, endpos))
        newarr.append(list(newseq))

    return (np.array(newarr), r)



def drawMiniAlignment(arr, log, nams, outfile, typ, dpi, title, width, height,
                      markup=False, markupdict=None):
    ali_height, ali_width = np.shape(arr)

    if  typ == 'nt':
        cD = getNtColours()
    elif typ == 'aa':
        cD = getAAColours()

    f = plt.figure(figsize=(width, height), dpi=dpi)

    a = f.add_subplot('111')
    a.set_xlim(0, ali_width)
    a.set_ylim(0, ali_height)
    xpoints = range(0, ali_width)
    j = 0.5
    for row in arr:
        cols = []
        for s in row:
            cols.append(cD[s])
        a.scatter(xpoints, [j]*ali_width, color=cols, s=0.75**2, marker="|", zorder=2)
        j += 1

    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)
    a.yaxis.set_ticks_position('none')
    a.yaxis.set_visible(False)

    if title:
        f.suptitle(title)
    for t in a.get_xticklabels():
        t.set_fontsize(8)

    if markup:
        if "remove_short" in markupdict:
            colour = "#d7ddf2"
            for nam in markupdict['remove_short']:
                i = nams.index(nam)
                a.hlines((i - 0.5), 0, ali_width, color=colour, lw=0.75, zorder=0)
        if "remove_insertions" in markupdict:
            colour = "#fece88"
            a.vlines(list(markupdict['remove_insertions']), 0, ali_height, color=colour, lw=0.75, zorder=1)
    f.savefig(outfile, dpi=dpi)


def main():
    parser = argparse.ArgumentParser(
            description='''Improve a multiple sequence alignment''')

    parser.add_argument("--infile", dest='infile', type=str,
                        help='path to input alignment')
    parser.add_argument("--outfile_stem", dest='outfile_stem', type=str,
                        help="stem for output files (including path)")
    parser.add_argument("--insertion_min_size", dest="insertion_min_size",
                        type=int, default=3,
                        help="miniumum size insertion to remove")
    parser.add_argument("--insertion_max_size", dest="insertion_max_size",
                        type=int, default=300,
                        help="miniumum size insertion to remove")
    parser.add_argument("--insertion_min_flank", dest="insertion_min_flank",
                        type=int, default=5,
                        help="minimum number of bases on either side of deleted insertions")
    parser.add_argument("--remove_min_length", dest="remove_min_length",
                        type=int, default=50)
    parser.add_argument("--crop_ends_mingap", dest='crop_ends_mingap',
                        type=int, default=10,
                        help="minimum gap size to crop from ends")
    parser.add_argument("--dpi", dest="plot_dpi",
                        type=int, default=300,
                        help="dpi for plots")
    parser.add_argument("--plot_format", dest="plot_format",
                        type=str, default='png',
                        help="plot format (png or svg)")
    parser.add_argument("--plot_width", dest="plot_width",
                        type=int, default=5,
                        help="width for plots (inches)")
    parser.add_argument("--plot_height", dest="plot_height",
                        type=int, default=3,
                        help="height for plots (inches)")
    parser.add_argument("--remove_insertions", dest="remove_insertions",
                        action="store_true")
    parser.add_argument("--remove_short", dest="remove_short",
                        action="store_true")
    parser.add_argument("--remove_gaponly", dest="remove_gaponly",
                        action="store_false")
    parser.add_argument("--crop_ends", dest="crop_ends",
                        action="store_true")
    parser.add_argument("--plot_input", dest="plot_input",
                        action="store_true")
    parser.add_argument("--plot_output", dest="plot_output",
                        action="store_true")
    parser.add_argument("--plot_markup", dest="plot_markup",
                        action="store_true")

    args = parser.parse_args()

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    handler = logging.FileHandler("log.log")
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    log.addHandler(handler)

    log.info("Initial parameters: %s" % str(args))

    # convert the input fasta file into an array and make a list of
    # sequence names so the order can be maintained
    arr, nams = FastaToArray(args.infile)

    print (arr)
    # store a copy of the original array
    orig_arr = copy.copy(arr)
    orig_nams = copy.copy(nams)

    # make a dictionary to store the changes made
    markupdict = dict()

    removed_seqs = set()
    removed_cols = set()

    # detect if the sequence is amino acids or nucleotides
    typ = seqType(arr)

    if typ == 'aa':
        log.info("Amino acid alignment detected")
    else:
        log.info("Nucleotide alignment detected")

    if args.remove_insertions:
        arr, r = removeInsertions(arr,
                                  log,
                                  args.insertion_min_size,
                                  args.insertion_max_size,
                                  args.insertion_min_flank)
        markupdict['remove_insertions'] = r
        removed_cols = removed_cols | r

    if args.remove_short:
        arr, r = removeTooShort(arr, log, args.remove_min_length, nams)
        markupdict['remove_short'] = r
        removed_seqs = removed_seqs | r

    if args.crop_ends:
        arr, r = cropEnds(arr, log, nams, args.crop_ends_mingap)
        markupdict['crop_ends'] = r

    if args.remove_gaponly:
        arr, r = removeGapOnly(arr, log)
        markupdict['remove_gaponly'] = r

    if args.plot_input:
        outf = "%s_input.%s" % (args.outfile_stem, args.plot_format)
        drawMiniAlignment(orig_arr, log, nams, outf, typ, args.plot_dpi,
                          args.outfile_stem, args.plot_width, args.plot_height)

    if args.plot_output:
        outf = "%s_output.%s" % (args.outfile_stem, args.plot_format)
        drawMiniAlignment(arr, log, nams, outf, typ, args.plot_dpi,
                          args.outfile_stem, args.plot_width, args.plot_height)

    if args.plot_markup:
        outf = "%s_markup.%s" % (args.outfile_stem, args.plot_format)
        drawMiniAlignment(orig_arr, log, nams, outf, typ, args.plot_dpi,
                          args.outfile_stem, args.plot_width, args.plot_height,
                          markup=True, markupdict=markupdict)



    outfile = "%s_parsed.fasta" % (args.outfile_stem)
    writeOutfile(outfile, arr, nams)


if __name__ == "__main__":
    main()
