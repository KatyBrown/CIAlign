#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
import matplotlib.font_manager
import sys
try:
    import CIAlign.palettes as palettes
except ImportError:
    import palettes
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


def FastaToArray(infile, log=None, outfile_stem=None):
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
                continue  # todo: test!
            if line[0] == ">":
                seqs.append([s.upper() for s in seq])
                nams.append(nam)
                seq = []
                nam = line.replace(">", "")
            else:
                if len(nams) == 0:
                    if log:
                        log.error(formatErrorMessage)
                    print(formatErrorMessage)
                    exit()
                seq += list(line)
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam)
    arr = np.array(seqs[1:])
    return (arr, nams[1:])


def getPalette(palette='CBS'):
    '''
    Generates a dictionary which assigns a name to each colour using a colour
    blindness safe palette, generated using
    https://medialab.github.io/iwanthue/
    Parameters
    ----------
    palette: str
        The ID of the palette to be used, currently only colour blind safe
        (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are names of colours and
        values are hexadecimal codes for colours
    '''
    if palette.lower() == 'cbs':
        p = palettes.CBSafe()

    return (p)


def getAAColours(pal='CBS'):
    '''
    Generates a dictionary which assigns a colour to each amino acid.
    Based on the "RasmMol" amino colour scheme and the table here:
        http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
    approximated using a CB safe palette generated using
    https://medialab.github.io/iwanthue/

    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter amino acid codes and
        values are hexadecimal codes for colours
    '''
    pal = getPalette(palette=pal)
    return {'D': pal['red_aa'],
            'E': pal['red_aa'],
            'C': pal['yellow_aa'],
            'M': pal['yellow_aa'],
            'K': pal['blue_aa'],
            'R': pal['blue_aa'],
            'S': pal['orange_aa'],
            'T': pal['orange_aa'],
            'F': pal['midblue_aa'],
            'Y': pal['midblue_aa'],
            'N': pal['cyan_aa'],
            'Q': pal['cyan_aa'],
            'G': pal['lightgrey_aa'],
            'L': pal['green_aa'],
            'V': pal['green_aa'],
            'I': pal['green_aa'],
            'A': pal['darkgrey_aa'],
            'W': pal['purple_aa'],
            'H': pal['paleblue_aa'],
            'P': pal['peach_aa'],
            'X': pal['black'],
            '-': pal['white'],
            'B': pal['tan_aa'],
            'Z': pal['tan_aa'],
            'J': pal['tan_aa'],
            '*': pal['white'],
            'U': pal['tan_aa'],
            'O': pal['tan_aa']
            }


def getNtColours(pal='CBS'):
    '''
    Generates a dictionary which assigns a colour to each nucleotide (plus grey
    for "N" and white for "-")
    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter nucleotide codes and
        values are hexadecimal codes for colours
    '''
    pal = getPalette(palette=pal)
    return {'A': pal['green_nt'],
            'G': pal['yellow_nt'],
            'T': pal['red_nt'],
            'C': pal['blue_nt'],
            'N': pal['grey_nt'],
            "-": pal['white'],
            "U": pal['red_nt'],
            "R": pal['grey_nt'],
            "Y": pal['grey_nt'],
            "S": pal['grey_nt'],
            "W": pal['grey_nt'],
            "K": pal['grey_nt'],
            "M": pal['grey_nt'],
            "B": pal['grey_nt'],
            "D": pal['grey_nt'],
            "H": pal['grey_nt'],
            "V": pal['grey_nt'],
            "X": pal['grey_nt']}


def getMarkupColours(pal='CBS'):
    '''
    Generates a dictionary which assigns a colour to each markup type
    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are CIAlign cleaning function names and
        values are hexadecimal codes for colours
    '''
    pal = getPalette(palette=pal)
    return {'remove_insertions': pal['remove_insertions'],
            'crop_ends': pal['crop_ends'],
            'remove_gaponly': pal['remove_gaponly'],
            'remove_short': pal['remove_short'],
            'remove_divergent': pal['remove_divergent']}


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
    plat = sys.platform
    if plat == "win32" or plat == "cygwin":
        raise RuntimeError("The list fonts function is currently \
                            unavailable in Windows")
    with warnings.catch_warnings():
        # Don't raise warnings for missing glyphs
        warnings.filterwarnings('ignore', message='Glyph')
        try:
            matplotlib.font_manager._rebuild()
        except AttributeError:
            pass
        flist = matplotlib.font_manager.get_fontconfig_fonts()
        flist2 = set()
        checkglyphs = [108, 112, 65, 71, 84, 67]
        for fname in flist:
            try:
                F = matplotlib.font_manager.FontProperties(fname=fname)
                font = matplotlib.font_manager.get_font(fname)
                L = font.get_charmap()
                # this tests if matplotlib can actually render "ACTG" in this
                # font
                x = 0
                for glyph in checkglyphs:
                    if glyph in L:
                        x += 1
                if x == len(checkglyphs):
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
        i = 0
        for i, fname in enumerate(flist2):
            a.text(0.7, i, "ACTG", fontdict={'name': fname, 'size': 14})
            a.text(0, i, fname, fontsize=10)

        a.text(0.7, i+1, "Sample", fontsize=10, fontweight='bold')
        a.text(0, i+1, "Font Name", fontsize=10, fontweight='bold')

        a.set_axis_off()

        f.tight_layout()
        f.savefig(outfile, dpi=200, bbox_inches='tight')
