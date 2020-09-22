#! /usr/bin/env python
'''
This script contains functions used to calculate various alignment statistics
which were used to benchmark CIAlign.
'''
import configargparse
import numpy as np
import copy
import itertools
import sys
# temporary until I sort out my $PATH
sys.path.insert(0, "/home/katy/CIAlign/CIAlign")
import utilityFunctions
# put this back eventually
#try:
#    import CIAlign.utilityFunctions as utilityFunctions
#except ImportError:
#    import utilityFunctions


def alignment_to_matrix(arr):
    '''
    Converts an alignment matrix stored as a numpy array into a
    matrix of integers, representing
    the position of each nucleotide in the sequence.

    This matrix is used in the comparison of two alignments, to recognise
    either columns or parts of aligned residues which are common or different
    between the two.

    Gaps ("-") and residues removed with CIAlign ("!") are represented as
    zeroes.

    e.g

    ATCGG : 1 2 3 4 5
    A-GCC : 1 0 2 3 4

    A-TCGG : 1 0 2 3 4 5
    AG--CC : 1 2 0 0 3 4

    Parameters
    ----------
    arr: np.array
        Numpy array containing the alignment represented as a 2D matrix, where
        dimension 1 is sequences and dimension 2 is columns

    Returns
    -------
    pos: np.array
        Numpy array of integers with the same dimensions as the input array
        showing the position of each non-gap residue in the original alignment
        as described above.
    '''
    # Each position is the cumulative number of non-gap residues to the left
    # of this position
    pos = np.cumsum(arr != "-", 1)
    # replace gaps in the input with 0
    pos[arr == "-"] = 0
    # replace ! (CIAlign removed residues) in the input with 0
    pos[arr == "!"] = 0
    pos = pos.astype(int)
    return(pos)


def find_removed_cialign(removed_file, arr, nams):
    '''
    Reads the "_removed.txt" file generated by CIAlign to determine
    what CIAlign has removed from the original alignment.

    Replaces nucleotides removed by CIAlign with "!" in the array representing
    the alignment so that it is still possible to compare these alignments
    with uncleaned alignments in terms of knowing which columns and pairs
    of residues are aligned.

    ! characters are always counted as mismatches in comparisons between
    alignments.

    Parameters
    ----------
    removed_file: str
        Path to a CIAlign _removed.txt log file
    arr: np.array
        Numpy array containing the alignment represented as a 2D matrix, where
        dimension 1 is sequences and dimension 2 is columns
    nams: list
        List of names in the original alignment, in the same order as in the
        input and the sequence array (these should always be the same).

    Returns
    -------
    cleanarr:
        2D numpy array containing the alignment represented as a 2D matrix,
        where dimension 1 is sequences and dimension 2 is columns, with
        residues removed by CIAlign represented as !
        Fully removed sequences are removed from this array.
    cleannams:
        List of names in the output alignment, with any sequences fully
        removed by CIAlign removed.
    '''
    # Read the CIAlign _removed.txt log file
    lines = [line.strip().split("\t")
             for line in open(removed_file).readlines()]
    # Make an empty dictionary
    D = {x: set() for x in nams}
    for line in lines:
        func = line[0]
        ids = line[-1].split(",")
        # for crop_ends and remove_insertions columns are removed so keep
        # track of column numbers as integers
        if func in ['crop_ends', 'remove_insertions']:
            ids = [int(x) for x in ids]
        # crop_ends is only applied to some sequences so also
        # keep track of sequence names
        if func == "crop_ends":
            nam = line[1]
            D[nam] = D[nam] | set(ids)
        # no need to remove insertions from sequences which were removed
        # completely later
        elif func == "remove_insertions":
            for nam in nams:
                if D[nam] != "removed":
                    D[nam] = D[nam] | set(ids)
        # remove divergent and remove short remove the whole sequence
        elif func in ["remove_divergent", "remove_short"]:
            for nam in ids:
                D[nam] = "removed"

    # make copies of the arrays (because I'm never quite sure when
    # python makes links rather than copies)
    cleannams = copy.copy(nams)
    cleannams = np.array(cleannams)
    cleanarr = copy.copy(arr)

    # iterate through everything that has been changed
    for nam, val in D.items():
        which_nam = np.where(cleannams == nam)[0][0]
        # remove the removed sequences from the array
        if val == "removed":
            cleannams = np.append(cleannams[:which_nam],
                                  cleannams[which_nam + 1:])
            cleanarr = np.vstack([cleanarr[:which_nam],
                                  cleanarr[which_nam+1:]])
            # remove them from the input temporarily just to keep the shapes
            # the same
            arr = np.vstack([arr[:which_nam], arr[which_nam+1:]])
        else:
            # replace column substitutions with !
            which_pos = np.array(sorted(list(val)))
            if len(which_pos) != 0:
                cleanarr[which_nam, which_pos] = "!"
    # sometimes gaps are removed - make these gaps in the output rather than
    # !s
    cleanarr[arr == "-"] = "-"
    return (cleanarr, cleannams)


def get_POARS(arr, nams):
    '''
    Converts an alignment matrix from alignment_to_matrix into a set of
    strings representing "pairs of aligned residues" as described in
    Lassman and Sonnhammer 2005 (DOI: 10.1186/1471-2105-8-S5-S9).

    Each string represents a single pair of residues and the sequences they
    originate from, delimited by _.
    Pairs containing gaps ("-") and residues removed
    by CIAlign (!) are not considered.

    e.g.
    >X
    ATC : 1 2 3
    >Y
    A-G : 1 0 2
    >Z
    AT- : 1 2 0

    POARS are X1-Y1, X1-Z1, X2-Z2, X3-Y2
    (strings would be X_Y_1_1, X_Z_1_1, X_Z_2_2, X_Y_3_2)

    Parameters
    ----------
    arr: np.array
        Numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue at this position in the sequence.
    nams: list
        List of sequence names in the same order as the rows of the sequence
        arry

    Returns
    -------
    POARS: set
        Set of strings representing all the pairs of aligned residues in the
        array as w_x_y_z where w is the name of sequence 1, x the name of
        sequence 2, y the index of the residue in sequence 1, z the
        index of the residue in sequence 2.
    '''
    nrows, ncols = np.shape(arr)
    POARS = set()
    # take every rowise combination of sequences in the alignment
    for i, j in itertools.combinations(np.arange(nrows), 2):
        # find the sequence names
        nam1 = nams[i]
        nam2 = nams[j]
        # take every column of this pair of sequences
        for k in np.arange(ncols):
            # get the indicies of the residues at these positions
            POAR = list(arr[[i, j], k])
            # if they are not gaps and haven't been removed with CIAlign
            if 0 not in POAR:
                # store the POAR
                POARS.add("%s_%s_%s_%s" % (nam1, nam2, POAR[0], POAR[1]))
    return (POARS)


def sum_of_pairs_overlap_score(POARS_benchmark, POARS_test):
    '''
    Calculate the "sum of pairs" score to compare two alignments described in
    Thompson et al. 1999 (DOI: 10.1093/nar/27.13.2682) and the "overlap score"
    described in Lassman and Sonnhammer 2005 (doi.org/10.1093/nar/gki1020)

    The sum of pairs score is calculated as the number of POARs present in
    both alignments divided by the number of POARs in the reference alignment.

    The overlap score is simular but is the number of POARS present in
    both alignments divided by the mean number of POARS in the two alignmnents.

    Sum of pairs score requires on alignment to be the reference and shows
    how "correct" a test alignment is.

    Overlap score can compare any two alignments and is a measure of how
    similar they are.

    Parameters
    ----------
    POARS_benchmark: set
        Set of strings representing all the pairs of aligned residues in the
        first alignment as w_x_y_z where w is the name of sequence 1,
        x the name of sequence 2, y the index of the residue in sequence 1,
        z the index of the residue in sequence 2.
    POARS_test: set
        Set of strings representing all the pairs of aligned residues in the
        first alignment as w_x_y_z where w is the name of sequence 1,
        x the name of sequence 2, y the index of the residue in sequence 1,
        z the index of the residue in sequence 2.


    Returns
    -------
    sum_of_pairs_score: float
        The sum of pairs score comparing the test alignment with the benchmark
        alignment
    overlap_score: float
        The overlap score comparing the two alignments
    '''

    # Find POARs present in both alignments
    POARS_both = POARS_benchmark & POARS_test

    # sum of pairs = N POARS in both / N POARS in benchmaark
    sum_of_pairs_score = len(POARS_both) / len(POARS_benchmark)

    # overlap = N POARS in both / mean POARs across both
    mean_POARS = (len(POARS_benchmark) + len(POARS_test)) / 2
    overlap_score = len(POARS_both) / mean_POARS

    return (sum_of_pairs_score, overlap_score)


def column_score(pos_benchmark, pos_test):
    '''
    Calculates the column score described in Thompson et al. 1999
    (DOI: 10.1093/nar/27.13.2682) to
    compare a test alignment with a benchmark "true" alignment.

    The score is the number of columns which are present in both
    alignments divided by the total number of columns in the true alignment.

    A correct column can be at a different position in the alignment
    but must consist of the residues from the same position as in the original
    alignment

    Based on a matrix showing the position of each nucleotide in the alignment
    generated using the alignment_to_matrix function.

    e.g.
    benchmark alignment
    ATCGG : 1 2 3 4 5
    A-GCC : 1 0 2 3 4

    test alignment
    A-TCGG : 1 0 2 3 4 5
    AG--CC : 1 2 0 0 3 4

    1-1, 4-3 and 5-4 are in both - so 3 columns are present in both alignments.

    The benchmark alignment has 5 total columns so the score would be 3/5

    Where alignments have been cleaned with CIAlign, ! symbols are inserted to
    allow this indexing to be accurate

    e.g.
    benchmark alignment
    ATCGG : 1 2 3 4 5
    A-GCC : 1 0 2 3 4

    cialign cleaned test alignment
    !-T-CGG: 0 0 2 0 3 4 5
    !G---CC: 0 2 0 0 0 3 4

    All zero columns are excluded

    2-0, 4-3 and 5-4 are in both so the score would be 3/5

    Parameters
    ----------
    pos_benchmark: np.array
        Numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue at this position in the sequence
        for the benchmark alignment (or the first for overlap score)
    pos_test: np.array
        Numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue at this position in the sequence
        for the test alignment (or the second for overlap score)

    Returns
    -------
    column_score: float
        The column score for the alignment - the number of columns present
        in both the benchmark and the test divided by the number of columns
        in the benchmark

    '''
    # don't count columns which are all zeroes
    pos_benchmark = pos_benchmark[:, pos_benchmark.sum(0) != 0]
    pos_test = pos_test[:, pos_test.sum(0) != 0]

    # make a set of strings where each string represents a column in the
    # benchmark
    S_orig = set(np.apply_along_axis(
        lambda x: "_".join(x), 0, pos_benchmark.astype(str)))

    # make a set of strings where each string represents a column in the
    # test
    S_new = set(np.apply_along_axis(
        lambda x: "_".join(x), 0, pos_test.astype(str)))

    # calculate the column score
    column_score = len(S_new & S_orig) / np.shape(pos_benchmark)[1]
    return (column_score)


def format_alignment(ali, cleaned=False, cialign_removed=None):
    '''
    Converts the alignment in the path ali to a numpy array of integers
    showing the cumulative number of non-gap residues prior to the residue
    at this position in the sequence, with characters removed by CIAlign
    excluded.
    Runs the FastaToArray function from utilityFunctions, converts to upper
    case, runs find_removed_cialign and alignment_to_matrix.

    Parameters
    ----------
    ali: str
        path to multiple sequence alignment in FASTA format. If the alignment
        has been cleaned with CIAlign this should be the CIAlign input, not
        the output
    cleaned: bool
        True if the alignment has been cleaned with CIAlign, otherwise False
    cialign_removed: str
        path to CIAlign _removed.txt file for the alignment

    Returns
    -------
    arr: np.array
        Numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue
        at this position in the sequence, with characters removed by CIAlign
        excluded
    nams: list
        List of sequence names in the same order as the rows of the sequence
        array.
    '''
    # Convert alignment into arrays
    arr, nams = utilityFunctions.FastaToArray(ali)

    # make everything upper case so this doesn't affect the score
    arr = np.char.upper(arr)

    # if the alignment has been cleaned with CIAlign, update the array
    # to contain !s for positions which have been removed
    if cleaned:
        arr, nams = find_removed_cialign(cialign_removed, arr, nams)

    arr = alignment_to_matrix(arr)

    return (arr, nams)


def calculate_alignment_scores(arr_1, arr_2, nams_1, nams_2):
    '''
    Calculates alignment scores above to compare a pair of alignments.
    Where a reference alignment is required, alignment_1 is assumed to be the
    reference.

    Parameters
    ----------
    arr_1: np.array
        2D numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue at this position in the sequence,
        with characters removed by CIAlign excluded, for the benchmark
        or first alignment.
    arr_2: np.array
        2D numpy array of integers showing the cumulative number of non-gap
        residues prior to the residue at this position in the sequence,
        with characters removed by CIAlign excluded, for the test or second
        alignment.
    nams_1: np.array
        List of sequence names in the same order as the rows of the sequence
        array for the benchmark or first alignment
    nams_2: np.array
        List of sequence names in the same order as the rows of the sequence
        array for the test or second alignment

    Returns
    -------
    sum_of_pairs_score: float
        The sum of pairs score comparing the test alignment with the benchmark
        alignment
    overlap_score: float
        The overlap score comparing the two alignments

    column_score: float
        The column score for the alignment - the number of columns present
        in both the benchmark and the test divided by the number of columns
        in the benchmark
    '''
    POARS_1 = get_POARS(arr_1, nams_1)
    POARS_2 = get_POARS(arr_2, nams_2)
    this_column_score = column_score(arr_1, arr_2)
    this_sum_of_pairs, this_overlap = sum_of_pairs_overlap_score(POARS_1,
                                                                 POARS_2)
    return (this_column_score, this_sum_of_pairs, this_overlap)


def main():
    parser = configargparse.ArgumentParser(
                description='Calculate scores to compare two alignments',
                add_help=False)

    parser.add("--alignment_1", dest='alignment_1', type=str,
               help='Path to first or benchmark alignment')

    parser.add("--alignment_2", dest='alignment_2', type=str,
               help='Path to second or test alignment')

    parser.add("--cialign_1", dest="cialign_1", type=str,
               default=None,
               help="Path to CIAlign _removed.txt for first alignment")
    parser.add("--cialign_2", dest="cialign_2", type=str,
               default=None,
               help="Path to CIAlign _removed.txt for second alignment")
    parser.add("--outfile", dest="outfile", type=str,
               default="scores.out",
               help="Path to output file")

    args = parser.parse_args()
    if args.cialign_1 is None:
        arr_1, nams_1 = format_alignment(args.alignment_1)
    else:
        arr_1, nams_1 = format_alignment(args.alignment_1,
                                         True, args.cialign_1)
    if args.cialign_2 is None:
        arr_2, nams_2 = format_alignment(args.alignment_2)
    else:
        arr_2, nams_2 = format_alignment(args.alignment_2,
                                         True, args.cialign_2)
    scores = calculate_alignment_scores(arr_1, arr_2, nams_1, nams_2)
    print(scores)


if __name__ == "__main__":
    main()
