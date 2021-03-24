#!/usr/bin/env python3
import numpy as np
import matplotlib
try:
    import CIAlign.cropSeq as cropSeq
except ImportError:
    import cropSeq
matplotlib.use('Agg')


def cropEnds(arr, nams, relativePositions, rmfile, log, mingap, redefine_perc):
    '''
    Removes poorly aligned ends from a multiple sequence alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object
    mingap: int
        Minimum gap size to crop from ends.

    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    r: dict
        A dictionary with sequence names as keys and tuples as values,
        where tuple[0] is a list of positions which have been removed at the
        beginning of the sequence and tuple[1] is a list of positions which
        have been removed at the end of the sequence
    '''
    out = open(rmfile, "a")
    newarr = []
    r = dict()
    for i, row in enumerate(arr):
        start, end = cropSeq.determineStartEnd(row, nams[i], log, mingap, redefine_perc)
        start = max(start - 1, 0)
        end = end + 1
        newseq = "-" * start + "".join(
                row[start:end]) + "-" * (len(row) - end)
        newseq = np.array(list(newseq))
        s = sum(newseq != row)

        if s != 0:
            nam = nams[i]
            non_gap_start = sum(row[0:start] != "-")
            non_gap_end = sum(row[end:] != "-")
            # list of positions between 0 and start which are not gaps
            startpos = np.where(row[0:start] != "-")[0]
            # list of positions from end to end of the sequences
            # which are not gaps
            endpos = np.where(row[end:] != "-")[0] + end

            rel_startpos = np.array(relativePositions)[startpos]
            rel_endpos = np.array(relativePositions)[endpos]

            if non_gap_start != 0:
                log.info("Removed %i bases from start of %s" % (
                        non_gap_start, nam))
                startpos_str = [str(x) for x in rel_startpos]
                out.write("crop_ends\t%s\t%s\n" % (nam,
                                                   ",".join(startpos_str)))
            if non_gap_end != 0:
                log.info("Removed %i bases from end of %s" % (
                        non_gap_end, nam))
                endpos_str = [str(x) for x in rel_endpos]
                out.write("crop_ends\t%s\t%s\n" % (nam,
                                                   ",".join(endpos_str)))

            r[nam] = ((rel_startpos, rel_endpos))
            # r[nam] = ((startpos, endpos))
        newarr.append(list(newseq))
    out.close()
    return (np.array(newarr), r)


def removeDivergent(arr, nams, rmfile, log, percidentity=0.75):
    '''
    Remove sequences which don't have the most common non-gap residue at
    > percidentity non-gap positions

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object
    percidentity: float
        Minimum percentage identity to majority to not be removed

    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    r: set
        A set of names of sequences which have been removed
    '''
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')
    j = 0
    keep = []
    for a in arr:
        i = 0
        y = 0
        t = 0
        if sum(a != "-") != 0:
            # for non-gap positions only (in this sequence)
            for base in a:
                if base != "-":
                    others = arr[:, i]
                    # only look at non-gap positions in the other sequences
                    others = others[others != "-"]
                    counts = np.unique(others, return_counts=True)
                    # find the most common
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
    newarr = arr[keep, :]
    r = set(np.array(nams)[np.invert(keep)])
    if len(r) != 0:
        log.info("Removing divergent sequences %s" % (", ".join(list(r))))
        out.write("remove_divergent\t%s\n" % (",".join(list(r))))
    out.close()
    return (newarr, r)


def removeInsertions(arr, relativePositions, rmfile, log,
                     min_size, max_size, min_flank):
    '''
    Removes insertions which are not present in the majority of sequences.
    Insertions are removed if they are between min_size and
    max_size residues and more than half of sequences with >= min_flank
    residues on either side of the insertion have the insertion.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment.
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object


    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    r: set
        A set of column numbers of sequences which have been removed
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment, minus
        the columns removed using this function.
    '''
    log.info("Removing insertions\n")
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')
    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column
    sums = sum(boolarr)
    # run a sliding window along the alignment and check for regions
    # which have higher coverage at the ends of the window than in the
    # middle - store these in put_indels
    put_indels = []
    for size in range(min_size, max_size, 1):
        for i in range(0, len(sums)+1 - size, 1):
            these_sums = sums[i:i+size]
            # take the number of non gap positions
            # for each column in this window
            ns = np.arange(i, i+size)
            left = these_sums[0]
            right = these_sums[-1]
            # record sites with lower coverage than their flanking seqs
            # and less than mincov total coverage
            x = set(ns[(these_sums < left) &
                       (these_sums < right)])
            if len(x) >= min_size:
                put_indels += list(x)
    put_indels = set(put_indels)
    put_indels = np.array(sorted(list(put_indels)))
    # for the putative indels, check if there are more sequences
    # with a gap at this position (but with sequence on either side)
    # than with no gap (but with sequence on either side)
    rmpos = set()
    absolutePositions = set()
    i = 0

    for p in put_indels:
        left = boolarr[:, :p]
        right = boolarr[:, p+1:]
        pcol = boolarr[:, p]
        pcol_nongaps = pcol
        pcol_gaps = np.invert(pcol)
        leftsum = np.sum(left, 1)
        rightsum = np.sum(right, 1)
        covers_region = (sum((pcol_nongaps) & (
                leftsum >= min_flank) & (rightsum >= min_flank)))
        lacks_region = (sum((pcol_gaps) & (
                leftsum >= min_flank) & (rightsum >= min_flank)))
        if lacks_region > covers_region:
            absolutePositions.add(p)
        i += 1

    # make a list of positions to remove
    rm_relative = set()
    for n in absolutePositions:
        rm_relative.add(relativePositions[n])
    for n in rm_relative:
        relativePositions.remove(n)
    # for n in absolutePositions:
    #     relativePositions.remove(n)
    rmpos = np.array(list(absolutePositions))

    keeppos = np.arange(0, len(sums))
    keeppos = np.invert(np.in1d(keeppos, rmpos))
    if len(rmpos) != 0:
        rmpos_str = [str(x) for x in rm_relative]
        log.info("Removing sites %s" % (", ".join(rmpos_str)))
        out.write("remove_insertions\t%s\n" % (",".join(rmpos_str)))
    out.close()
    arr = arr[:, keeppos]
    # return (arr, set(rmpos), relativePositions)
    return (arr, set(rm_relative), relativePositions)


def removeTooShort(arr, nams, rmfile, log, min_length):
    '''
    Removes sequences (rows) with fewer than min_length non-gap positions from
    the alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object
    min_length: int
        Minimum length sequence to keep.

    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    rmnames: set
         A set of names of sequences which have been removed

    '''
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')
    if len(arr) != 0:
        arrT = arr.transpose()
        sums = sum(arrT != "-")
        arr = arr[sums > min_length]
        rmnames = set(np.array(nams)[sums <= min_length])
        if len(rmnames) != 0:
            log.info("Removing too short sequences %s" % (
                ", ".join(list(rmnames))))
            out.write("remove_too_short\t%s\n" % ",".join(list(rmnames)))
    else:
        rmnames = set()
    out.close()
    return (arr, rmnames)


def removeGapOnly(arr, relativePositions, rmfile, log):
    '''
    Removes gap only columns from the alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment.
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object

    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    r: set
        A set of column numbers of sequences which have been removed
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment, minus
        the columns removed using this function.
    '''
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')
    if len(arr) != 0:
        sums = sum(arr == "-")
        absolutePositions = set(np.where(sums == len(arr[:, 0]))[0])
        rmpos = []
        for n in absolutePositions:
            rmpos.append(relativePositions[n])
        # remove deleted columns from relativePositions
        for n in rmpos:
            relativePositions.remove(n)
        rmpos = set(rmpos)
        arr = arr[:, sums != len(arr[:, 0])]
        if len(rmpos) != 0:
            rmpos_str = [str(x) for x in rmpos]
            log.info("Removing gap only sites %s" % (
                    ", ".join(rmpos_str)))
            out.write("remove_gap_only\t%s\n" % (",".join(rmpos_str)))
    else:
        rmpos = set()

    out.close()
    return (arr, rmpos, relativePositions)
