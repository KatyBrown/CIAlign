#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import cropSeq


def cropEnds(arr, nams, rmfile, log, mingap):
    out = open(rmfile, "a")
    newarr = []
    r = dict()
    for i, row in enumerate(arr):
        start, end = cropSeq.determineStartEnd(row, mingap)
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
                out.write("Removed %i bases from start of %s" % (non_gap_start, nam))
                out.write('\n')
            if non_gap_end != 0:
                log.info("Removed %i bases from end of %s" % (non_gap_end, nam))
                out.write("Removed %i bases from end of %s" % (non_gap_end, nam))
                out.write('\n')
            startpos = np.where(row[0:start] != "-")[0]
            endpos = np.where(row[end:] != "-")[0] + end
            r[nam] = ((startpos, endpos))
        newarr.append(list(newseq))
    out.close()
    return (np.array(newarr), r)


def removeBadlyAligned(arr, nams, rmfile, log, percidentity=0.9):
    '''
    Remove sequences which don't have the most common non-gap residue at
    > percidentity non-gap positions
    '''
    j = 0
    keep = []
    for a in arr:
        i = 0
        y = 0
        t = 0
        if sum(a != "-") != 0:
            for base in a:
                if base != "-":
                    others = arr[:, i]
                    others = others[others != "-"]
                    counts = np.unique(others, return_counts=True)
                    mc = counts[0][counts[1] == max(counts[1])][0]
                    if base == mc:
                        y += 1
                    t += 1
                i += 1
            if y / t > percidentity:
                keep.append(True)
            else:
                keep.append(False)
        else:
            keep.append(False)
        j += 1
    keep = np.array(keep)
    newarr = arr[keep,:]
    r = set(np.array(nams)[np.invert(keep)])
    return (newarr, r)


def removeInsertions(arr, relativePositions, rmfile, log,
                     min_size, max_size, min_flank):
    '''
    Removes insertions of size between min_size and
    max_size which have lower coverage than their flanking reg ions, if at least
    twice as many other sequences have the flanking regions
    but not the insertions.
    '''
    log.info("Removing insertions\n")
    out = open(rmfile, "a")
    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column
    sums = sum(boolarr)
    # run a sliding window along the alignment and check for regions
    # which have higher coverage at the ends of the window than in the
    # middle - store these in put_indels
    put_indels = set()
    for size in range(min_size, max_size, 1):
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
    # for the putative indels, check if there are more sequences
    # with a gap at this position (but with sequence on either side)
    # than with no gap (but with sequence on either side)
    rmpos = set()
    absolutePositions = set()
    for p in put_indels:
        left = arr[:,:p]
        right = arr[:,p+1:]
        pcol = arr[:,p]
        pcol_nongaps = pcol != "-"
        pcol_gaps = pcol == "-"
        leftsum = sum(left.T != "-")
        rightsum = sum(right.T != "-")
        covers_region = (sum((pcol_nongaps) & (leftsum >= min_flank) & (rightsum >= min_flank)))
        lacks_region = (sum((pcol_gaps) & (leftsum >= min_flank) & (rightsum >= min_flank)))
        if lacks_region > covers_region:
            absolutePositions.add(p)
    # make a list of positions to keep
    for n in absolutePositions:
        rmpos.add(relativePositions[n])
    for n in absolutePositions:
        relativePositions.remove(n)
    rmpos = np.array(list(rmpos))
    keeppos = np.arange(0, len(sums))
    #keeppos = np.invert(np.in1d(keeppos, rmpos)) # I think here we need the absolute positions? not sure though
    keeppos = np.invert(np.in1d(keeppos, rmpos))
    log.info("Removing sites %s" % (", ".join([str(x) for x in rmpos])))
    out.write("Removing sites %s" % (", ".join([str(x) for x in rmpos])))
    out.write('\n')
    out.close()
    arr = arr[:, keeppos]
    return (arr, set(rmpos))


def removeTooShort(arr, nams, rmfile, log, min_length):
    '''
    Removes sequences with fewer than min_length non-gap positions from
    the alignment.
    '''
    if len(arr) != 0:
        arrT = arr.transpose()
        sums = sum(arrT != "-")
        arr = arr[sums > min_length]
        rmnames = set(np.array(nams)[sums <= min_length])
        log.info("Removing sequences %s" % (", ".join(list(rmnames))))
    else:
        rmnames = set()
    return (arr, rmnames)


def removeGapOnly(arr, relativePositions, rmfile, log):
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')

    if len(arr) != 0:
        print (arr.shape)
        sums = sum(arr == "-")
        absolutePositions = set(np.where(sums == len(arr[:,0]))[0])
        rmpos = []
        for n in absolutePositions:
            rmpos.append(relativePositions[n])
        # remove deleted columns from relativePositions
        for n in rmpos:
            relativePositions.remove(n)
        rmpos = set(rmpos)
        arr = arr[:, sums != len(arr[:,0])]
        log.info("Removing gap only sites %s" % (", ".join([str(x) for x in rmpos])))
        out.write("Removing gap only sites %s" % (", ".join([str(x) for x in rmpos])))
        out.write('\n')
    else:
        rmpos = set()

    out.close()
    return (arr, rmpos)







