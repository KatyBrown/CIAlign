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
    Calculates a frequency matrix to use to generate a PWM, showing the
    background frequency of the base at each position.
    
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
        ff = np.repeat(1 / len(core), len(core))
    elif freqtype == 'calc':
        ff = np.sum(PFM.values, 1) / np.sum(PFM.values)
    elif freqtype == 'calc2':
        assert PFM2 is not None, "To use PWM frequency method calc2 a second \
PFM must be provided"
        ff = np.sum(PFM2.values, 1) / np.sum(PFM2.values)
    elif isinstance(freqtype, dict):
        for res in core:
            assert res in freqtype, "Residue %s not found in \
frequency dictionary when calculating PWM frequency" % res
        ff = [freqtype[res] for res in core]
    else:
        raise RuntimeError ("Frequency type not recognised for PWM frequency")
    
    ffmat = np.reshape(np.repeat(ff, np.shape(PFM)[1]), np.shape(PFM))
    return (ffmat)


def getAlpha(alphatype, PFM, freq):
    if alphatype == 'calc':
        alpha = np.array([f * np.sqrt(np.sum(PFM)) for f in freq])
    elif alphatype is None:
        alpha = np.full(np.shape(PFM), 0)
    elif isinstance(alphatype, int):
        alpha = np.full(np.shape(PFM), alphatype)
    return (alpha)


def makePFM(arr, typ):
    core = getCoreRes(typ)
    sums = [np.sum(arr == res, 0) for res in core]
    sum_mat = pd.DataFrame(sums, index=core)
    return (sum_mat.round(4))

def makePPM(PFM, alpha):
    PPM = np.divide(PFM + alpha, np.sum(PFM + alpha))
    return (PPM.round(4))

def makePWM(PPM, freq):
    PWM = np.log2(np.divide(PPM, freq))
    return (PWM.round(4))