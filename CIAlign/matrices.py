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


def getFreq(freqtype, log, typ, PFM, PFM2=None):
    '''
    Calculates a frequency matrix to use to generate a position probability
    matrix, showing the background frequency of the base at each position -
    the expected frequency in a random sequence.
    
    There are various ways to calculate this, 3 are currently implemented:
    
    If freqtype is 'equal' it is assumed that all residues are equally common,
    so the frequency of each will be 0.25 for nucleotides and 0.05 for
    amino acids.
    
    If freqtype is "calc" the frequency is calculated using the PFM, by
    dividing the total number of occurances of each residue by the total
    number of residues.
    
    If freqtype is "calc2" a second PFM is used to calculate the frequency,
    for example this could be a full alignment if the PWM you are calculating
    is only for part of the alignment.

    Parameters
    ----------
    freqtype: str
        Should be "equal", "calc" or "calc2" as
        described above.
    log: logging.Logger
        An open log file object
    typ: str
        nt for nucleotide or aa for amino acid
    PFM: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in the alignment, calculated with the makePFM function.
    PFM2: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in a second alignment which will be used for the background
        frequency calculation if freqtype == 'calc2'.
    
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
        log.info(
            "Assuming PPM background frequencies are equal at %.4f" % ff[0])
    elif freqtype == 'calc':
        # Number of each residue / total number of residues
        ff = np.sum(PFM.values, 1) / np.sum(PFM.values)
        log.info(
            "Calculating PPM background frequencies based on alignment")
    elif freqtype == 'calc2':
        # Check a second PFM exists
        assert PFM2 is not None, "To use PWM frequency method calc2 a second \
PFM must be provided"
        # As above but for the second PFM
        ff = np.sum(PFM2.values, 1) / np.sum(PFM2.values)
        log.info(
            "Calculating PPM background frequencies based on \
             secondary alignment")
    else:
        raise RuntimeError ("Frequency type not recognised for PWM frequency")
    
    # The calculations above are per nucleotide, repeat for each column
    ffmat = np.reshape(np.repeat(ff, np.shape(PFM)[1]), np.shape(PFM))
    return (ffmat)


def getAlpha(alphatype, log, PFM, freq, alphaval=1.0):
    '''
    Calculates the alpha parameter used as a pseudocount when creating a
    position probability matrix to avoid zero values.
    There are various ways to calculate this, currently two are implemented.
    
    If alphatype is "calc", alpha is calculated as
    frequency(base) * (square root(number of rows in alignment)), as described
    in Dave Tang's blog here:
    https://davetang.org/muse/2013/10/01/position-weight-matrix/, which
    recreates the method used in doi.org/10.1038/nrg1315
    
    If alphatype is "user" the value is provided by the user as alphaval.

    Parameters
    ----------
    alphatype: str
        "calc", to calculate alpha as described above or "user" to provide a
        value of alpha as alphaval
    log: logging.Logger
        An open log file object
    PFM: np.array
        Array showing the frequency of each nucleotide at each position
        (column) in the alignment, calculated with the makePFM function.
    freq: np.array
        Array showing the background frequency of the residue at each position,
        calculated using the getFreq function.
    alphaval: float
        User defined value of alpha  

    Returns
    -------
    alpha: np.array
        An array with the same dimensions as the PFM, showing the alpha value
        for each residue.
    '''

    if alphatype == 'calc':
        # Calculate as discussed in the blog post linked above
        alpha = np.array([f * np.sqrt(np.sum(PFM)) for f in freq])
        log.info("Calculating alpha value based on background frequencies")
    elif isinstance(alphatype, int) or isinstance(alphatype, float):
        # If an integer is provided, use the integer
        # Reshape to match the PFM
        alpha = np.full(np.shape(PFM), alphaval)
        log.info("Using user defined alpha value of %.4f" % alphaval)
    else:
        raise RuntimeError ("Frequency type not recognised for PWM frequency")
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
    of each residue in each column of the alignment, normalised by background
    frequency and including pseudocounts.
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
