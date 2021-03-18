#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


def replaceUbyT(arr):
    '''
    Replaces all Us by Ts in the alignment.

    Parameters
    ----------
    arr: np.array
        2D numpu array with the alignment.

    Returns
    -------
    arr: np.array
        2D numpy array of sequences with Ts instead of Us.
    '''

    arr = np.where(arr == "U", "T", arr)
    return (arr)


def unAlign(arr):
    '''
    Removes all gaps from the alignment.

    Parameters
    ----------
    arr: np.array
        2D numpu array with the alignment.

    Returns
    -------
    arr: np.array
        2D numpy array of sequences without any gaps.
    '''

    arr = np.where(arr == "-", "", arr)
    return (arr)


def FastaToArray(infile, log, outfile_stem=None):
    '''
    Convert an alignment into a numpy array.

    Parameters
    ----------
    infile: string
        path to input alignment file in FASTA format
    log: logging.Logger
        An open log file object

    Returns
    -------
    arr: np.array
        2D numpy array in the same order as fasta_dict where each row
        represents a single column in the alignment and each column a
        single sequence.
    nams: list
        List of sequence names in the same order as in the input file
    '''

    formatErrorMessage = "The MSA file needs to be in FASTA format."
    nams = []
    seqs = []
    nam = ""
    seq = ""
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue # todo: test!
            if line[0] == ">":
                seqs.append([s.upper() for s in seq])
                nams.append(nam)
                seq = []
                nam = line.replace(">", "")
            else:
                if len(nams) == 0:
                    log.error(formatErrorMessage)
                    print(formatErrorMessage)
                    exit()
                seq += list(line)
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam)
    arr = np.array(seqs[1:])
    return (arr, nams[1:])


def getAAColours():
    '''
    Generates a dictionary which assigns a colour to each amino acid.
    Based on the "RasmMol" amino colour scheme and the table here:
        http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
    (plus grey for "X" and white for "-")

    Parameters
    ----------
    None

    Returns
    -------
    dict
        Dictionary where keys are single letter amino acid codes and
        values are hexadecimal codes for colours
    '''
    return {'D': '#E60A0A',
            'E': '#E60A0A',
            'C': '#E6E600',
            'M': '#E6E600',
            'K': '#145AFF',
            'R': '#145AFF',
            'S': '#FA9600',
            'T': '#FA9600',
            'F': '#3232AA',
            'Y': '#3232AA',
            'N': '#00DCDC',
            'Q': '#00DCDC',
            'G': '#EBEBEB',
            'L': '#0F820F',
            'V': '#0F820F',
            'I': '#0F820F',
            'A': '#C8C8C8',
            'W': '#B45AB4',
            'H': '#8282D2',
            'P': '#DC9682',
            'X': '#b2b2b2',
            '-': '#FFFFFF00',
            'B': '#b2b2b2',
            'Z': '#b2b2b2',
            'J': '#b2b2b2',
            '*': '#FFFFFF00'
            }


def getNtColours():
    '''
    Generates a dictionary which assigns a colour to each nucleotide (plus grey
    for "N" and white for "-")
    Parameters
    ----------
    None

    Returns
    -------
    dict
        Dictionary where keys are single letter nucleotide codes and
        values are hexadecimal codes for colours
    '''
    return {'A': '#1ed30f',
            'G': '#f4d931',
            'T': '#f43131',
            'C': '#315af4',
            'N': '#b2b2b2',
            "-": '#FFFFFF',
            "U": '#f43131',
            "R": '#b2b2b2',
            "Y": '#b2b2b2',
            "S": '#b2b2b2',
            "W": '#b2b2b2',
            "K": '#b2b2b2',
            "M": '#b2b2b2',
            "B": '#b2b2b2',
            "D": '#b2b2b2',
            "H": '#b2b2b2',
            "V": '#b2b2b2',
            "X": '#b2b2b2'}


def writeOutfile(outfile, arr, nams, removed, rmfile=None):
    '''
    Writes an alignment stored in a numpy array into a FASTA file.

    Parameters
    ----------
    outfile: str
        Path to FASTA file where the output should be stored
    arr: np.array
        Numpy array containing the cleaned alignment
    nams: list
        List of nams of sequences in the input alignment
    removed: set
        Set of names of sequences which have been removed
    rmfile: str
        Path to file used to log sequences and columns which have been removed

    Returns
    -------
    None

    '''
    out = open(outfile, "w")
    i = 0
    for nam in nams:
        if nam not in removed:
            out.write(">%s\n%s\n" % (nam, "".join(list(arr[i]))))
            i += 1
    out.close()


def seqType(arr):
    '''
    Detects if an alignment is of nucleotides or amino acids using pre-built
    dictionarys of amino acid and nucleotide codes.
    Checks if arr contains characters that are not in the dictionary (not
    IUPAC)

    Parameters
    ----------
    arr: np.array
        Numpy array containing the alignment

    Returns
    -------
    str
    'aa' for amino acid and 'nt for nucleotide
    '''
    nt_count = 0
    aa_count = 0
    for seq in arr:
        nucs = set(list(getNtColours().keys()))
        aas = set(list(getAAColours().keys()))
        n = 0
        a = 0
        x = 0
        for s in seq:
            s = s.upper()
            if s in nucs:
                n += 1
            if s in aas:
                a += 1
            if s not in aas and s not in nucs:
                x += 1
        ch = 0
        if n == len(seq):
            nt_count += 1
            ch += 1
        if a == len(seq):
            aa_count += 1
            ch += 1
        if ch == 0:
            print("Unknown nucleotides or amino acids detected.\
                  Please fix your MSA.")
            exit()

    if nt_count == len(arr):
        return "nt"
    if aa_count == len(arr):
        return "aa"
    print("MSA type couldn't be established. Please fix your MSA.")
    exit()


def updateNams(nams, removed_seqs):
    '''
    Takes nams, a list of sequence names in the input file, and removed_seqs,
    a set of sequences which have been removed from the file and subtracts
    the removed sequences from nams while maintaining the order.

    Parameters
    ----------
    nams: list
        A list of sequence names in the input file
    removed_seqs: set
        Set of sequence names to remove

    Returns
    -------
    nams2: list
        A list of sequence names which have not been removed
    '''
    nams2 = []
    for nam in nams:
        if nam not in removed_seqs:
            nams2.append(nam)
    return (nams2)


def checkArrLength(arr, log):
    '''
    Checks the shape of the array containing the alignment to ensure that it
    all the sequences haven't been removed by parsing and that all the
    sequences have the same number of columns

    Parameters
    -----------
    arr: np.array
        Numpy array containing the multiple sequence alignment
    log: logging.Logger
        An open log file object
    Returns
    -------
    None
    '''
    emptyAlignmentMessage = """Error: Parsing your alignment with these \
settings has removed all of the sequences."""
    differentLengthsMessage = """Error: The sequences in your alignment are \
not all the same length."""
    if 0 in np.shape(arr):
        log.error(emptyAlignmentMessage)
        print(emptyAlignmentMessage)
        exit()
    if len(np.shape(arr)) == 1:
        log.error(differentLengthsMessage)
        print(differentLengthsMessage)
        exit()


def listFonts(outfile):
    '''
    Generates an image file containing a sample of all the fonts matplotlib
    has access to on your system, which can be used to choose a font for a
    text-based sequence logo.
    Not parameterised - the same image will be generated every time.

    Parameters
    ----------
    outfile: str
        Path to the file in which to save the image

    Returns
    -------
    None
    '''
    matplotlib.font_manager._rebuild()
    flist = matplotlib.font_manager.get_fontconfig_fonts()
    flist2 = set()
    for fname in flist:
        try:
            F = matplotlib.font_manager.FontProperties(fname=fname)
            font = matplotlib.font_manager.get_font(fname)
            L = font.get_charmap()
            # this tests if matplotlib can actually render "ACTG" in this
            # font
            if 108 in L:
                g = F.get_name()
                flist2.add(g)
        except RuntimeError:
            # Some of the fonts seem not to have a name? Ignore these.
            pass

    flist2 = sorted(list(flist2))[::-1]
    f = plt.figure(figsize=(5, len(flist2) / 4), dpi=200)
    a = f.add_subplot(111)
    a.set_ylim(0, len(flist2))
    a.set_xlim(0, 1)
    a.text(-0.1, -1, "*Fonts shown as [] cannot be displayed with CIAlign")
    for i, fname in enumerate(flist2):
        a.text(0.7, i, "ACTG", fontdict={'name': fname, 'size': 14})
        a.text(0, i, fname, fontsize=10)

    a.text(0.7, i+1, "Sample", fontsize=10, fontweight='bold')
    a.text(0, i+1, "Font Name", fontsize=10, fontweight='bold')
    a.set_axis_off()
    f.tight_layout()
    f.savefig(outfile, dpi=200, bbox_inches='tight')
