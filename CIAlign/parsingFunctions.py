#!/usr/bin/env python3
import numpy as np
import matplotlib
try:
    import CIAlign.cropSeq as cropSeq
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.insertions as insertions
except ImportError:
    import cropSeq
    import utilityFunctions
    import insertions
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

    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column
    sums = sum(boolarr)
    
    height, width = np.shape(arr)

    low_coverage = insertions.findLowCoverage(boolarr,
                                              sums, height, width,
                                              min_size, max_size)

    put_indels = insertions.getPutativeIndels(boolarr, sums, width,
                                              low_coverage, min_size,
                                              max_size)
    good_indels = insertions.findGoodInsertions(put_indels,
                                                boolarr,
                                                min_size, max_size, min_flank)
    absolutePositions = insertions.finalCheck(good_indels, min_size, max_size)

    arr, relativePositions, rm_relative = utilityFunctions.removeColumns(
        absolutePositions, relativePositions, arr, log, rmfile,
        "remove_insertions")
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

    if len(arr) != 0:
        sums = sum(arr == "-")
        absolutePositions = set(np.where(sums == len(arr[:, 0]))[0])
        
        arr, relativePositions, rmpos = utilityFunctions.removeColumns(
            absolutePositions, relativePositions, arr, log, rmfile,
            "remove_gaponly")
    else:
        rmpos = set()

    return (arr, rmpos, relativePositions)
