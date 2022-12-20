#!/usr/bin/env python3
import numpy as np
import math
import copy


def findLowCoverage(boolarr, sums, height, width, min_size, max_size,
                    min_flank, min_perc):
    '''
    Finds regions of the alignment with coverage <50%, as insertions gaps must
    be found in the majority of sequences to be removed -
    regions with >=50% coverage can never meet this criteria.
    Generates a list of arrays, each of which is the indices of a group of
    consecutive positions with coverage in <50% of columns. The groups
    must have length >= min_size and < max_size.

    Parameters
    ----------
    boolarr: np.array
        Boolean numpy array the same shape as the alignment where residues are
        True and gaps are False.
    sums: np.array
        Integer array showing the total number of non-gaps in each column
        of the alignment
    height: int
        Alignment height (number of sequences)
    width: int
        Alignment width (number of columns)
    min_size: int
        Minimum size insertions to remove
    max_size: int
        Maximum size insertions to remove
    min_perc: float
        Remove insertions which are present in less than this
        proportion of sequences. Default: 0.5
    Returns
    -------
    puts: np.array
        Integer array of arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences.
    '''

    # Set the minimum threshold - the column shouldn't have coverage in
    # >thresh of rows
    thresh = math.floor(height) * min_perc

    # Create a boolean array with True if columns have >thresh coverage
    # and False otherwise
    low_cov = (sums <= thresh)

    # Calculate the cumulative sum of the number of columns with <50% coverage
    # - the numbers only increment for regions with coverage below the
    # threshold
    cumarr = np.cumsum(low_cov)

    # Find the positions where the cumulative array changes from 0 to 1 - the
    # breakpoints between coverage above and below the threshold
    puts = np.where(np.diff(cumarr) == 1)[0] + 1
    diffs = np.where(np.diff(puts) != 1)[0] + 1
    splits = np.split(puts, diffs)
    # Add a buffer few residues to low coverage regions - for some complex
    # insertions the iterative process
    # later somehow clips off the last residue or two
    splits2 = []
    for split in splits:
        if len(split) != 0:
            s = copy.deepcopy(split)
            for i in np.arange(s[-1], s[-1]+6):
                s = np.append(s, i)
            splits2.append(s)
    return (np.array(splits2, dtype='object'))


def getPutativeIndels(boolarr, sums, width, puts,
                      min_size, max_size, min_flank, min_perc):
    '''
    Narrows down the low coverage regions to keep only low coverage regions
    flanked by higher coverage regions. Regions which don't meet these
    criteria are split into smaller regions if this then allows them
    to meet the criteria.

    Parameters
    ----------
    boolarr: np.array
        Boolean numpy array the same shape as the alignment where residues are
        True and gaps are False.
    sums: np.array
        Integer array showing the total number of non-gaps in each column
        of the alignment
    width: int
        Alignment width (number of columns)
    puts: np.array
        Integer array of arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences.
    min_size: int
        Minimum size insertions to remove
    max_size: int
        Maximum size insertions to remove
    min_flank: int
        Minimum number of bases on either side of an insertion
        to classify it as an insertion.
    min_perc: float
        Remove insertions which are present in less than this
        proportion of sequences. Default: 0.5
    Returns
    -------
    pp: list
        List of numpy integer arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences, flanked by columns with higher coverage.

    '''
    # Store the results
    pp = []
    ufD = dict()
    width = np.shape(boolarr)[1]
    leftsum = np.cumsum(boolarr, 1)
    rightsum = (np.cumsum(boolarr[:, ::-1], 1)[:, ::-1])
    # Iterate through the list of low coverage regions
    for put in puts:
        if len(put) != 0:
            # Variable to track the change in the start position
            # The maximum size of the deletion is the size of the current
            # low coverage region
            # The deletion can't be bigger than the alignment
            ms = min([len(put), (width-put[0]-2), (max_size - 3)])
            if ms < 0:
                break
            # Start with the biggest
            current_size = ms
            # The end of the insertion to test initially is the end of the
            # current
            # low coverage region

            # Iterate until you get to the end of the LC region or the current
            # size is 0
            while (current_size >= min_size):
                current_start = put[0]
                current_end = current_start + current_size
                while (current_end <= put[-1]):
                    # Take the region inside the putative deletion
                    wi = sums[current_start:current_end]
                    # Take the indicies immediately on either side
                    left_pos = current_start - 1
                    right_pos = current_end

                    if right_pos < width:
                        # Get the number of gaps in the columns on either side
                        before_pos = np.sum(boolarr[:, left_pos])
                        after_pos = np.sum(boolarr[:, right_pos])
                        # If everything within the put deletion has lower
                        # coverage than either side

                        if np.all(wi < before_pos) & np.all(wi < after_pos):
                            good = []
                            bad = []
                            # Iterate through the individual positions
                            for p in np.arange(current_start, current_end):

                                # Take the current pos
                                pcol_nongaps = boolarr[:, p]
                                pcol_gaps = np.invert(pcol_nongaps)

                                # Count the non-gap positions on either side
                                # of the put deletion in each row
                                ls = leftsum[:, p-1]
                                rs = rightsum[:, p+1]

                                # Check that more columns have the insertion
                                # than don't (from those with enough flanking
                                # residues)
                                covers_region = (sum((pcol_nongaps) & (
                                        ls >= min_flank) & (rs >= min_flank)))
                                lacks_region = (sum((pcol_gaps) & (
                                        ls >= min_flank) & (rs >= min_flank)))

                                # Check that the majority of rows lack
                                # the insertion
                                if lacks_region + covers_region > 0:
                                    prop_with_insertion = covers_region / (
                                        lacks_region + covers_region)

                                    if prop_with_insertion <= min_perc:
                                        good.append(p)
                                    else:
                                        bad.append(p)
                                else:
                                    bad.append(p)

                                # Save the flank info for the next iteration
                                fl = leftsum[:, p-1]
                                fr = rightsum[:, p+1]
                                uf = fl, fr
                                ufD[p] = uf
                            # If all the positions are OK
                            if len(bad) == 0:
                                # Save this deletion
                                pp.append(np.arange(current_start,
                                                    current_end+1))
                                # move on to the next one
                                current_start += current_size
                                current_end = current_start + current_size
                            else:
                                # try the next bigger size
                                current_start += 1
                                current_end = current_start + current_size
                        else:
                            # try the next smaller size
                            current_start += 1
                            current_end = current_start + current_size
                    else:
                        current_start += 1
                        current_end = current_start + current_size
                else:
                    # try the next smaller size
                    current_start += 1
                    current_end = current_start + current_size

                # Reduce the current size
                current_size -= 1
    return (pp)


def finalCheck(pp, sums, min_size, max_size):
    '''
    Check all the insertions are between min_size and max_size and flanked
    by a higher coverage region and remove any which don't meet these criteria.

    Parameters
    ----------
    absolutePositions : set
        Set of positions which to potentially remove as insertions
    min_size : min_size
        Minimum size insertions to remove
    max_size : max_size
        Maximum size insertions to remove

    Returns
    -------
    keep : set
        Final set of positions to remove
    '''
    absolutePositions = set()
    for put in pp:
        for p in put[:-1]:
            absolutePositions.add(p)
    abso = np.array(sorted(list(absolutePositions)))
    # Find the consecutive positions and split the array at these points
    diffs = np.where(np.diff(abso) != 1)[0] + 1
    splits = np.array(np.split(abso, diffs), dtype=object)
    # Filter to keep only the insertions which are within the size limit
    keep = splits[np.array([(len(r) >= min_size) &
                            (len(r) < max_size) for r in splits])]
    absolutePositions_final = set()
    for split in keep:
        current_start = split[0]
        current_end = split[-1]+1
        # Check that all positions within the insertion have higher coverage
        # than all positions immediately on either side
        xx = np.arange(current_start-1, current_end+1)
        cent = sums[xx[1:-1]]
        bef = sums[xx[0]]
        aft = sums[xx[-1]]
        if (np.all(cent < bef)) & (np.all(cent < aft)):
            absolutePositions_final = absolutePositions_final | set(split)
    return (absolutePositions_final)
