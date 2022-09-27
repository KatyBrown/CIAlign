#!/usr/bin/env python3
import numpy as np
import pandas as pd


def getCoreRes(typ):
    '''
    Returns the "standard" (non-ambiguous, non-gap) set of residues for
    either nucleotides or amino acids.
    
    Parameters
    ----------
    typ: str
        nt - nucleotide or aa - amino acid
    '''
    if typ == "nt":
        core = ['A', 'C', 'T', 'G']
    elif typ == "aa":
        core = ['D', 'E', 'C', 'M', 'K',
                'R', 'S', 'T', 'F', 'Y',
                'N', 'Q', 'G', 'L', 'V',
                'I', 'A', 'W', 'H', 'P']
    return (sorted(core))


def getFreq(freqtype, typ, PFM, PFM2=None):
    '''
    Calculates a frequency matrix to use to generate a position probability
    matrix, showing the background frequency of the base at each position -
    the expected frequency in a random sequence.
    
    There are various ways to calculate this, four are currently implemented:
    
    If freqtype is 'equal' it is assumed that all residues are equally common,
    so the frequency of each will be 0.25 for nucleotides and 0.05 for
    amino acids.
    
    If freqtype is "calc" the frequency is calculated using the PFM, by
    dividing the total number of occurances of each residue by the total
    number of residues.
    
    If freqtype is "calc2" a second PFM is used to calculate the frequency,
    for example this could be a full alignment if the PWM you are calculating
    is only for part of the alignment.
    
    If freqtype is a dictionary, the values of the dictionary are used
    as user defined frequencies for each residue.
    
    
    Parameters
    ----------
    freqtype: str or dict
        If freqtype is a string is should be "equal", "calc" or "calc2", as
        described above. If it is a dictionary the keys should be all possible
        residues and the values their frequencies as floats.
    typ: str
        nt for nucleotide or aa for amino acid
    PFM: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in the alignment, calculated with the makePFM function.
    PFM2: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in a second alignment which will be used for the background
        frequency calculation if freqtype == 'calc2'
    
    Returns
    -------
    freq: np.array
        An array with the same dimensions as the PFM showing the residue
        frequency at each position.
    '''
    core = getCoreRes(typ)
    if freqtype == 'equal':
        # 1 / number of possible residues
        ff = np.repeat(1 / len(core), len(core))
    elif freqtype == 'calc':
        # Number of each residue / total number of residues
        ff = np.sum(PFM.values, 1) / np.sum(PFM.values)
    elif freqtype == 'calc2':
        # Check a second PFM exists
        assert PFM2 is not None, "To use PWM frequency method calc2 a second \
PFM must be provided"

        # As above but for the second PFM
        ff = np.sum(PFM2.values, 1) / np.sum(PFM2.values)
    elif isinstance(freqtype, dict):
        # freqtype can also be a user provided dictionary
        for res in core:
            # Check all the possible residues have a frequency
            assert res in freqtype, "Residue %s not found in \
frequency dictionary when calculating PWM frequency" % res

        # Convert to a list in the right order
        ff = [freqtype[res] for res in core]
    else:
        raise RuntimeError ("Frequency type not recognised for PWM frequency")
    
    # The calculations above are per nucleotide, repeat for each column
    ffmat = np.reshape(np.repeat(ff, np.shape(PFM)[1]), np.shape(PFM))
    return (ffmat)


def getAlpha(alphatype, PFM, freq):
    '''
    Calculates the alpha parameter used as a pseudocount when creating a
    position probability matrix to avoid zero values.
    There are various ways to calculate this, currently two are implemented.
    
    If alphatype is "calc", alpha is calculated as
    frequency(base) * (square root(number of rows in alignment)), as described
    in Dave Tang's blog here:
    https://davetang.org/muse/2013/10/01/position-weight-matrix/, which
    recreates the method used in doi.org/10.1038/nrg1315
    
    If alphatype is None no alpha value is applied (alpha = 0)
    
    If alphatype is an integer or float, this value of alpha is used for
    all positions.


    Parameters
    ----------
    alphatype: str or int or None
        If alpha is a string it should be "calc", to calculate alpha as
        described above. If it is None alpha will be zero, if it is an integer
        or float this will be used as the value of alpha.
    PFM: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in the alignment, calculated with the makePFM function.
    freq: np.array
        Array showing the background frequency of the residue at each position,
        calculated using the getFreq function.
    
    Returns
    -------
    alpha: np.array
        An array with the same dimensions as the PFM, showing the alpha value
        for each residue.
    '''

    if alphatype == 'calc':
        # Calculate as discussed in the blog post linked above
        alpha = np.array([f * np.sqrt(np.sum(PFM)) for f in freq])
    elif alphatype is None:
        # If alphatype is None, just use 0
        # Reshape to match the PFM
        alpha = np.full(np.shape(PFM), 0)
    elif isinstance(alphatype, int) or isinstance(alphatype, float):
        # If an integer is provided, use the integer
        # Reshape to match the PFM
        alpha = np.full(np.shape(PFM), alphatype)
    return (alpha)


def makePFM(arr, typ):
    '''
    Make a position frequency matrix - a matrix showing the number of
    each residue in each column of the alignment. Gaps and ambiguous codes
    are ignored as most downstream software doesn't accept these.
    
    PFMs are described in:
        https://en.wikipedia.org/wiki/Position_weight_matrix
        10.1186/s12859-020-3348-6
        https://davetang.org/muse/2013/10/01/position-weight-matrix/
    This implementation reproduces all three of these examples exactly.
    
    Parameters
    ----------
    arr: np.array
        The alignment the PFM is for, stored as a numpy array
    typ: str
        nt for nucleotide or aa for amino acid
        
    Returns
    -------
    PFM: pd.DataFrame
        A position frequency matrix as a pandas Data Frame, where columns
        are columns in the alignment (numbered as 0 to n columns) and
        rows are residues - either A, C, G T or the 20 standard amino acids
        in alphabetical order.   
    '''
    core = getCoreRes(typ)
    sums = [np.sum(arr == res, 0) for res in core]
    sum_mat = pd.DataFrame(sums, index=core)
    return (sum_mat.round(4))

def makePPM(PFM, alpha):
    '''
    Make a position probability matrix - a matrix showing the frequency
    of each residue in each column of the alignment, normalised by alignment
    size. 
    Alpha is a pseudocount - a small number added to avoid zero counts. It
    can be provided directly as an integer or calculated as a transformation
    of the input data as described in the getFreq function docstring.
    
    Gaps and ambiguous codes
    are ignored as most downstream software doesn't accept these.
    
    PFMs are described in:
        https://en.wikipedia.org/wiki/Position_weight_matrix
        10.1186/s12859-020-3348-6
    This implementation reproduces these examples.
    
    Parameters
    ----------
    arr: np.array
        The alignment the PFM is for, stored as a numpy array
    typ: str
        nt for nucleotide or aa for amino acid
        
    Returns
    -------
    PPM: pd.DataFrame
        A position probability matrix as a pandas Data Frame, where columns
        are columns in the alignment (numbered as 0 to n columns) and
        rows are residues - either A, C, G T or the 20 standard amino acids
        in alphabetical order.
    '''
    PPM = np.divide(PFM + alpha, np.sum(PFM + alpha))
    return (PPM.round(4))


def makePWM(PPM, freq):
    '''
    Make a position weight matrix - a matrix which shows the log-likelihood
    ratio of observing character i at position j in a site
    compared with a random sequence (from 10.1186/s12859-020-3348-6).
    
    Freq is an array with the same dimensions as the PFM showing the expected
    residue frequency at each position in a random sequence, calculated using
    the getFreq function.
    
    Gaps and ambiguous codes
    are ignored as most downstream software doesn't accept these.
    
    PWMs are described in:
        https://en.wikipedia.org/wiki/Position_weight_matrix
        10.1186/s12859-020-3348-6
        https://davetang.org/muse/2013/10/01/position-weight-matrix/
    This implementation reproduces these examples.

    Parameters
    ----------
    arr: np.array
        The alignment the PFM is for, stored as a numpy array
    freq: np.array
        An array with the same dimensions as the PFM showing the expected
        residue frequency at each position in a random sequence, calculated
        using the getFreq function.
        
    Returns
    -------
    PWM: pd.DataFrame
        A position weight matrix as a pandas Data Frame, where columns
        are columns in the alignment (numbered as 0 to n columns) and
        rows are residues - either A, C, G T or the 20 standard amino acids
        in alphabetical order.     
    '''
    PWM = np.log2(np.divide(PPM, freq))
    return (PWM.round(4))