#!/usr/bin/env python3

import copy
import numpy as np


def cropDivergentPos(arr, min_prop_ident, min_prop_nongap, buffer,
                     start=True):
    '''
    Find the new start (if start=True) or end (if end=True) postion if
    the alignment is cropped to remove divergent flanking sequences.
    Find the index of the leftmost column in the first series of
    consecutive columns of length buffer where the proportion of
    non-gap residues doesn't fall below min_prop_nongap and the
    proportion of non-gap residues which are identical doesn't fall below
    min_prop_ident

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array
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
    int
        The index of the new start or end column for the alignment
    '''

    # Reverse the whole alignment (temporarily) to crop from the 3' end
    # Use copies of arrays for safety
    if start:
        this_arr = copy.copy(arr)
    else:
        this_arr = copy.copy(arr[:, ::-1])

    # Find the number of rows represented by the min_prop proportion
    min_count_nongap = int(np.shape(this_arr)[0] * min_prop_nongap)
    passfail = []
    # Iterate through the columns in the array, keep track of column index
    for i, col in enumerate(this_arr.T):
        # Isolate the non-gap rows in this column
        col_nongap = col[col != "-"]

        # Check if there are enough non-gap residues
        pass_nongap = len(col_nongap) >= min_count_nongap

        # If there are not enough non-gap residues this column is no
        # good
        if not pass_nongap:
            passfail.append(False)
            continue
        else:
            # Count how many of each nt/aa there are
            C = np.unique(col_nongap, return_counts=True)

            # Find which nt/aa is the most common
            which_max = C[1] == max(C[1])

            # Find how many times the most common nt/aa occurs
            max_col_count = C[1][which_max][0]
            max_col_prop = max_col_count / len(col_nongap)

            # Check if there are enough identical residues
            pass_ident = max_col_prop >= min_prop_ident

            # If there are not enough identical residues this column
            # is no good
            if not pass_ident:
                passfail.append(False)
                continue
            else:
                passfail.append(True)

            # Are all of the previous buffer columns True?
            buffer_passfail = passfail[-buffer:]

            # Once we meet the criteria we can stop
            if sum(buffer_passfail) == buffer:
                break

    # Find the position at the beginning of the buffer
    # Reverse if cropping the 3' end of the alignment
    if start:
        return (i - buffer + 1)
    else:
        return (np.shape(arr)[1] - i + buffer - 1)
