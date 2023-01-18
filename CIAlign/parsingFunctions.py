#!/usr/bin/env python3
import numpy as np
import matplotlib
try:
    import CIAlign.cropSeq as cropSeq
    import CIAlign.cropDiv as cropDiv
    import CIAlign.insertions as insertions
    import CIAlign.utilityFunctions as utilityFunctions
except ImportError:
    import cropSeq
    import cropDiv
    import insertions
    import utilityFunctions
matplotlib.use('Agg')


def cropEnds(arr, nams, relativePositions, rmfile, log, keeps,
             mingap, redefine_perc):
    '''
    Removes poorly aligned ends from a multiple sequence alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
    nams: list
        The names of the sequences in the alignment
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment.
    rmfile: str
        Path to a file in which to store a list of removed sequences
    log: logging.Logger
        An open log file object
    keeps: dict
        A dictionary listing sequences not to process for each function,
        where the keys are function names.
    mingap: int
        Minimum gap size to crop from ends.
    redefine_perc: float
        Proportion of the sequence length (excluding gaps) that is being
        checked for change in gap numbers to redefine start/end.

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
        start, end = cropSeq.determineStartEnd(row, nams[i], log,
                                               mingap, redefine_perc)
        # Find the sequences the user specified to keep in either
        # crop ends or all rowwise functions
        kk = set(keeps['crop_ends']) | set(keeps['all_rowwise'])
        # Don't redefine the start and end for the "retain" seqs
        if nams[i] in kk:
            start = 0
            end = len(row)
        else:
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


def removeDivergent(arr, nams, rmfile, log,
                    keeps, percidentity=0.75):
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
    keeps: dict
        A dictionary listing sequences not to process for each function,
        where the keys are function names.
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
    # Find which sequence names the user specified not to process
    # with this function
    kk = set(keeps['remove_divergent']) | set(keeps['all_rowwise'])

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
            # keep sequences in the kk array
            if y / t > percidentity or nams[j] in kk:
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
                     min_size, max_size, min_flank, min_perc):
    '''
    Removes insertions which are not present in > min_perc proportion of
    sequences.
    Insertions are removed if they are between min_size and
    max_size residues and more than min_perc proportion of sequences where
    >= min_flank residues on either side of the insertion have the insertion.

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
    min_size: int
        Only remove insertions >= this number of residues.
    max_size: int
        Only remove insertions <= this number of residues.
    min_flank: int
        Minimum number of bases on either side of an insertion
        to classify it as an insertion.
    min_perc: float
        Remove insertions which are present in less than this
        proportion of sequences. Default: 0.5


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
    parr = np.zeros([0, 0])
    i = 0
    rm_rel = set()
    absolutePositions = set()
    pufD = dict()
    abso = list(np.arange(np.shape(arr)[1]))

    # record which sites are not "-"
    boolarr = arr != "-"
    # array of the number of non-gap sites in each column
    sums = sum(boolarr)

    height, width = np.shape(arr)

    low_coverage = insertions.findLowCoverage(boolarr,
                                              sums, height, width,
                                              min_size, max_size,
                                              min_flank, min_perc)

    put_indels = insertions.getPutativeIndels(boolarr, sums, width,
                                              low_coverage, min_size,
                                              max_size, min_flank, min_perc)

    absolutePositions = insertions.finalCheck(put_indels, sums, min_size,
                                              max_size)

    arr, relativePositions, rm_relative = utilityFunctions.removeColumns(
        absolutePositions, relativePositions, arr, log, rmfile,
        "remove_insertions")

    return (arr, set(rm_relative), relativePositions)


def removeTooShort(arr, nams, rmfile, log, keeps, min_length):
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
    keeps: dict
        A dictionary listing sequences not to process for each function,
        where the keys are function names.
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
        # Find which sequence names the user specified not to process
        # with this function
        keeps_rs = np.in1d(nams, keeps['remove_short'])
        keeps_ar = np.in1d(nams, keeps['all_rowwise'])
        keeps_here = keeps_rs | keeps_ar
        # Artificially raise the number of non-gap residues in
        # the sequences listed in the "keeps" list so that they are
        # not removed but otherwise processed like all other sequences
        sums[keeps_here] += min_length
        arr = arr[(sums > min_length)]

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
            "remove_gap_only")
    else:
        rmpos = set()

    return (arr, rmpos, relativePositions)


def cropDivergent(arr, relativePositions, rmfile, log,
                  min_prop_ident, min_prop_nongap, buffer_size):
    '''
    Find the new start and end postion if the alignment is cropped
    to remove divergent flanking sequences.
    Find the index of the leftmost column in the first series of
    consecutive columns of length buffer where the proportion of
    non-gap residues doesn't fall below min_prop_nongap and the
    proportion of non-gap residues which are identical doesn't fall below
    min_prop_ident, repeat for right size of the alignment.

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
    min_prop_ident: float
        The minimum proportion of sequences which should have the same
        residue in each column
    min_prop_nongap: float
        The minimum proportion of sequences which should not be gaps in
        each column
    buffer: int
        The number of consecutive columns which should meet the min_prop_ident
        and min_prop_nongap criteria to pass filtering
    start: bool
        True - find the start position
        False - find the end position

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

    log.info("Cropping divergent sequence ends\n")
    out = open(rmfile, "a")
    if out.closed:
        print('file is closed')

    new_start = cropDiv.cropDivergentPos(arr,
                                         min_prop_ident,
                                         min_prop_nongap,
                                         buffer_size,
                                         start=True)
    new_end = cropDiv.cropDivergentPos(arr,
                                       min_prop_ident,
                                       min_prop_nongap,
                                       buffer_size,
                                       start=False)
    rmStart = np.arange(0, new_start)
    rmEnd = np.arange(new_end, np.shape(arr)[1])

    absolutePositions = np.append(rmStart, rmEnd)

    rm_relative = set()
    for n in absolutePositions:
        rm_relative.add(relativePositions[n])
    for n in rm_relative:
        relativePositions.remove(n)
    # for n in absolutePositions:
    #     relativePositions.remove(n)
    rmpos = np.array(list(absolutePositions))

    keeppos = np.arange(new_start, new_end)

    if len(rmpos) != 0:
        rmpos_str = [str(x) for x in rm_relative]
        log.info("Cropping divergent ends - keeping positions %i-%i" % (
            new_start, new_end))
        log.info("Removed positions: %s" % (",".join(rmpos_str)))
        out.write("crop_divergent\t%s\n" % (",".join(rmpos_str)))
    out.close()
    arr = arr[:, keeppos]
    # return (arr, set(rmpos), relativePositions)
    return (arr, set(rm_relative), relativePositions)
