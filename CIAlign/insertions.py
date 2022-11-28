#!/usr/bin/env python3
import numpy as np
import math


def findLowCoverage(boolarr, sums, height, width, min_size, max_size):
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
    
    Returns
    -------
    puts: np.array
        Integer array of arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences.
    '''
    # Set the minimum threshold - the column shouldn't have coverage in >50%
    # of rows
    thresh = math.floor(height/2)

    # Create a boolean array with True if columns have >50% coverage
    # and False otherwise
    low_cov = np.array(sums <= thresh)

    # Calculate the cumulative sum of the number of columns with <50% coverage
    # - the numbers only increment for regions with coverage below the
    # threshold
    cumarr = np.cumsum(low_cov)

    # Find the positions where the cumulative array changes from 0 to 1 - the
    # breakpoints between coverage above and below the threshold
    puts = np.where(np.diff(cumarr) == 1)[0] + 1

    diffs = np.where(np.diff(puts) != 1)[0] + 1
    splits = np.split(puts, diffs)
    return (splits)


def getPutativeIndels(boolarr, sums, width, puts,
                      min_size, max_size):
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

    Returns
    -------
    pp: list
        List of numpy integer arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences, flanked by columns with higher coverage.

    '''
    # Store the results
    pp = []
    for put in puts:
        current_start = put[0]
        current_size = min_size
        current_end = put[-1]
        # Test if the full length low coverage region is lower coverage than
        # the columns on either side
        # If not, keep making it bigger until you get a hit or reach max_size
        # - then move to the end of the hit and start again
        while (current_end <= put[-1]):
            while (current_size <= max_size):
                wi = sums[current_start:current_end]
                # Find the indices immediately either side of the low coverage
                # regions
                left_pos = current_start - 1
                right_pos = current_end + 1
                # If the low coverage region isn't at the end of the alignment
                if left_pos >= 0 and right_pos <= width-1:
    
                    # Get the number of gaps in the columns on either side
                    before = np.sum(boolarr[:, left_pos])
                    after = np.sum(boolarr[:, right_pos])
                    
                    if np.all(wi < before) & np.all(wi < after):
                        # if coverage is higher on either side of the region than
                        # inside it, it's already a candidate indel
                        pp.append(np.arange(current_start, current_end+2))
                        # move on to the next one
                        current_start += current_size
                        current_end = current_start + current_size
                        break
                    else:
                        # try the next bigger size
                        current_size += 1
                else:
                    # try the next bigger size
                    current_size += 1
                current_end = current_start + current_size
            
            current_start += 1
            # Reset the current size to be the maximum which will fit in
            # the remaining sequence
            # current_size =  np.min([len(put) - current_start - 2, max_size])
            current_size =  min_size

            # Reset the current end based on the current size and start
            current_end = current_start + current_size
    return (pp)


def findGoodInsertions(pp, boolarr, min_size, max_size, min_flank):
    '''
    Insertions which should be removed are those where the gap is present
    in more sequences than the inserted sequence in >50% of the total
    sequences spanning the insertion by a length of at least min_flank.
    
    Parameters
    ----------
    pp: list
        List of numpy integer arrays where each array contains the indices of
        a group of between min_size and max_size consecutive columns with
        gaps in >50% of sequences, flanked by columns with higher coverage.
    boolarr: np.array
        Boolean numpy array the same shape as the alignment where residues are
        True and gaps are False.
    min_size: int
        Minimum size insertions to remove
    max_size: int
        Maximum size insertions to remove
    min_flank: int
        Minimum number of bases on either side of an insertion to classify it
        as an insertion.

    Returns
    -------
    absolutePositions : set
        Set of positions to potentially remove

    '''
    absolutePositions = set()
    width = np.shape(boolarr)[1]
  #  print (width)
    for put in pp:
        pstart, pend = put[0], put[-1]
        left_lim = pstart - min_flank
        right_lim = pend + min_flank
        if left_lim >= 0 and right_lim <= width:
            left_seg = boolarr[:, left_lim:pstart]
            centre_seg = boolarr[:, pstart:pend]
            right_seg = boolarr[:, pend:right_lim]
            mm = np.max(np.sum(centre_seg, 0))
            leftsum = np.sum(left_seg, 0)
            rightsum = np.sum(right_seg, 0)
            if all(leftsum > mm) and all(rightsum > mm):
                for p in put[:-1]:
                    absolutePositions.add(p)
    return (absolutePositions)


def finalCheck(absolutePositions, min_size, max_size):
    '''
    Check all the insertions are between min_size and max_size and remove
    any which don't meet these criteria.

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
    abso = np.array(sorted(list(absolutePositions)))
    # Find the consecutive positions and split the array at these points
    diffs = np.where(np.diff(abso) != 1)[0] + 1
    splits = np.array(np.split(abso, diffs), dtype=object)
    # Filter to keep only the insertions which are within the size limit
    keep = splits[np.array([(len(r) >= min_size) &
                            (len(r) <= max_size) for r in splits])]
    # Flatten the array
    absolutePositions_final = set()
    for k in keep:
        absolutePositions_final = absolutePositions_final | set(k)
    return (absolutePositions_final)