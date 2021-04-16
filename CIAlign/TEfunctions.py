#! /usr/bin/env python
import copy
import numpy as np
import TEfunctions
try:
    import CIAlign.cropSeq as cropSeq
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.TEs as TEs
except ImportError:
    import cropSeq
    import utilityFunctions
    import TEs

def flankingMotifs():
    pass

def clusterTEs():
    pass



def cropDivergent(arr, relativePositions, rmfile, log,
                  cd_start, cd_end,
                  cd_min_prop_ident, cd_min_prop_nongap,
                  cd_buffer):
    '''
    Crops divergetn columns from either end of the alignment after
    identifying the first and/or last conserved region.
    
    Conserved regions are defined as a number (defined by the cd_buffer
    parameter of sequential columns which have a proportion of non-gap
    residues above a threshold (cd_min_prop_nongap) and
    a proportion of identical residues above a threshold (cd_min_prop_ident).

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
    cd_start: bool
        If true, crop divergent columns from the left side of the alignment
    cd_end: bool
        If true, crop divergent columns from the right side of the alignment
    cd_min_prop_ident: float
        The minimum proportion of sequences which should have the same
        residue in each column
    cd_min_prop_nongap: float
        The minimum proportion of sequences which should not be gaps in
        each column
    cd_buffer: int
        The number of consecutive columns which should meet the min_prop_ident
        and min_prop_nongap criteria to pass filtering
    
    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    rmpos: set
        A set of column numbers of sequences which have been removed
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment, minus
        the columns removed using this function.
    '''

    L = np.shape(arr)[1]
    if cd_start:
        new_start = TEs.cropDivergentPos(arr,
                                         min_prop_ident=cd_min_prop_ident,
                                         min_prop_nongap=cd_min_prop_nongap,
                                         buffer=cd_buffer,
                                         start=True)

    else:
        new_start = 0
    if cd_end:
        new_end = TEs.cropDivergentPos(arr,
                                       min_prop_ident=cd_min_prop_ident,
                                       min_prop_nongap=cd_min_prop_nongap,
                                       buffer=cd_buffer,
                                       start=False)
    else:
        new_end = L

    absolutePositions = list(np.arange(0, new_start))
    absolutePositions += list(np.arange(new_end, L))
    arr, relativePositions, rmpos = utilityFunctions.removeColumns(
            absolutePositions, relativePositions, arr, log, rmfile,
            "crop_divergent")
    return (arr, rmpos, relativePositions)
    
def cropDivergent_byrow():
    pass