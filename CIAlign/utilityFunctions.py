#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
import matplotlib.font_manager
import sys
import os
try:
    import CIAlign.palettes as palettes
except ImportError:
    import palettes
matplotlib.use('Agg')


def replaceUbyT(arr, rev):
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
    if not rev:
        arr = np.where(arr == "T", "U", arr)
    else:
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
                continue
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
                    exit(1)
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
    if palette.lower() == 'bright':
        p = palettes.Bright()
    return (p)


def getAAColours(palette='CBS'):
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
    pal = getPalette(palette=palette)
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


def getNtColours(palette='CBS'):
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
    pal = getPalette(palette=palette)
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
            'remove_gap_only': pal['remove_gap_only'],
            'remove_short': pal['remove_short'],
            'remove_divergent': pal['remove_divergent'],
            'crop_divergent': pal['crop_divergent'],
            'user': pal['user']}


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


def seqType(arr, log):
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
            log.error("Unknown nucleotides or amino acids detected.\
                       Please fix your MSA.")
            print("Unknown nucleotides or amino acids detected.\
                  Please fix your MSA.")
            exit(1)
    if nt_count == len(arr):
        return "nt"
    if aa_count == len(arr):
        return "aa"


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
        exit(1)
    if len(np.shape(arr)) == 1:
        log.error(differentLengthsMessage)
        print(differentLengthsMessage)
        exit(1)


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


def configRetainSeqs(retain, retainS, retainL, nams, fname, log, silent):
    '''
    Allows the user to specify sequences to keep regardless of whether
    they pass or fail the rowwise cleaning operation thresholds. This function
    works on the sequences specified for a single function or group of
    functions - either all rowwise functions, crop ends, remove divergent
    or remove short.

    Sequence names can be specified individually on the command line with
    --retain, --crop_ends_retain, --remove_divergent_retain,
    --remove_short_retain - these are listed in the retain variable for
    the specific function currently being processed.

    They can also be specified as a list in a text file with --retain_list,
    --crop_ends_retain_list etc, in which case the value is the path to
    the file.

    Finally they can be specified by searching each name for a character
    string, specified as --retain_str, --crop_ends_retain_str etc.

    retain, retainS and retainL are None if no sequences are specified.

    Parameters
    ----------
    retain: list
        List of sequence names to keep
    retainS: str
        Sequence names containing this string will be kept
    retainL: str
        Path to a text file containing a list of sequence names to keep
    nams: np.array
        Array containing the names of the sequences in the input fasta file
    fname: str
        The function currently being processed - all_rowwise, crop_ends,
        remove_divergent or remove_short, used for logging only here.
    log: logging.Logger
        An open log file object
    silent: bool
        True if CIAlign is run in silent mode - nothing will be printed to
        STDOUT

    Returns
    -------
    keeps: np.array
        An array containing all sequence names to be ignored by this
        function
    '''
    # First read the sequence names passed directly
    # the len check is for users with an inifile - can't specify None
    # so this will be an empty string
    if retain is not None and len(retain[0].strip()) != 0:
        keeps = set(retain)
    else:
        keeps = set()

    # If a file is specified, read the sequence names in the file
    if retainL is not None and len(retainL.strip()) != 0:
        # Raise an error if the file is not found
        if not os.path.exists(retainL):
            raise FileNotFoundError("""
List of sequences to retain %s not found""" % retainL)

        with open(retainL) as infile:
            for line in infile:
                keeps.add(line.strip())

    # If a string to match is specified
    if retainS is not None and len(retainS[0].strip()) != 0:
        # For each string
        for rs in retainS:
            rr = 0
            # check all the sequence names for this string
            for nam in nams:
                if rs in nam:
                    keeps.add(nam)
                    rr += 1
            if rr == 0:
                # Warn if there are no matches
                log.warning("""No sequence names matching "%s" were found"""
                            % rs)
                if not silent:
                    print("""Warning: No sequence names matching \
"%s" were found""" % rs)

    # Check all the specified sequence names were found, raise an error if
    # not
    if len(keeps & set(nams)) != len(keeps):
        raise RuntimeError("""
Some sequences listed to be retained were not found: %s""" % (" ".join(
                                                              keeps - set(
                                                                  nams))))

    # Convert the result to an array
    keeps_arr = np.array(sorted(list(keeps)))

    # Log the sequence names identified
    if len(keeps_arr) != 0:
        log.info("""The following sequences will not be processed with the \
%s function: %s""" % (fname, ", ".join(sorted(list(keeps)))))

    return (keeps_arr)


def updateStartEnd(start, end, removed):
    '''
    Updates the start and end positions of a subsection of an array taking
    into account positions that have been removed i.e. if the user
    wants columns 10:20 of the input FASTA but 5 and 15 have been removed,
    then in the updated alignment they would want 9:18.

    Parameters
    ----------
    start: int
        The start position in the unedited alignment
    end: int
        The end position in the unedited alignment
    removed: set
        Set of integers representing positions which have been removed from
        the alignment

    Returns
    -------
    newstart: int
        The updated start position
    newend: int
        The updated end position
    '''
    newstart = start
    newend = end
    if len(removed) != 0:
        for r in sorted(removed):
            if r < newstart:
                newstart -= 1
                newend -= 1
            elif r < newend:
                newend -= 1
    return (newstart, newend)


def removeColumns(rmAbsolute, relativePositions,
                  arr, log, rmfile, function_name, write=True):
    '''
    Function to remove a column from an array shared by removeInsertions,
    removeGapOnly. Removes the columns, calculates the
    positions of these columns in the input alignment, writes positions
    removed to the log file and the removed file.

    Parameters
    ----------
    rmAbsolute: list
        List of absolute positions (positions in current alignment) to
        be removed
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment.
    arr: np.array
        The alignment stored as a numpy array
    log: logging.Logger
        An open log file object
    rmfile: str
        Path to a file in which to store a list of removed sequences
    function_name: str
        The name of the function which removed these columns for logging
    Returns
    -------
    arr: np.array
        The cleaned alignment stored as a numpy array
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment, minus
        the columns removed using this function.
    r: set
        A set of column numbers of sequences which have been removed
    '''
    outrm = open(rmfile, "a")
    if outrm.closed:
        print('file is closed')
    # make a list of positions to remove
    rm_relative = set()
    for n in rmAbsolute:
        rm_relative.add(relativePositions[n])
    for n in rm_relative:
        relativePositions.remove(n)
    # for n in absolutePositions:
    #     relativePositions.remove(n)
    rmpos = np.array(list(rmAbsolute))

    keeppos = np.arange(0, np.shape(arr)[1])
    keeppos = np.invert(np.in1d(keeppos, rmpos))
    if len(rmpos) != 0 and write:
        rmpos_str = [str(x) for x in sorted(rm_relative)]
        log.info("Removing sites %s" % (", ".join(rmpos_str)))
        outrm.write("%s\t%s\n" % (function_name, ",".join(rmpos_str)))
    outrm.close()
    arr = arr[:, keeppos]
    return (arr, relativePositions, rm_relative)
