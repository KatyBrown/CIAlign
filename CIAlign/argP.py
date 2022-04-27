#!/usr/bin/env python3
import configargparse
import os.path
import numpy as np

try:
    import CIAlign.utilityFunctions as utilityFunctions
    from CIAlign._version import __version__
except ImportError:
    import utilityFunctions
    from _version import __version__


def float_range(mini, maxi):
    '''
    Defines a type for argparse of a float with a fixed range

    Parameters
    ----------
    mini: str or float
        A float (or a string which can be converted to a float) with the
        minumum valid value for the parameter.

    maxi: str or float
        A float (or a string which can be converted to a float) with the
        maximum valid value for the paramter.

    Returns
    -------
    float_range_checker: function
        A function to input as a type to argparse to check if the value
        provided is in the right range
    '''
    # Based on solution from @georg-w stack overflow qu. 55324449
    mini = float(mini)
    maxi = float(maxi)

    def float_range_checker(arg):
        try:
            f = float(arg)
        except ValueError:
            raise configargparse.ArgumentTypeError(
                "Must be a floating point number")
        if f < mini or f > maxi:
            raise configargparse.ArgumentTypeError(
                "Must be in range [%s .. %s]" % (mini, maxi))
        return (f)

    return (float_range_checker)


def int_range(mini, maxi, n_col):
    '''
    Defines a type for argparse of an integer with a fixed range, where
    the maximum value is a function of the number of columns

    Parameters
    ----------
    mini: str or int
        An integer (or a string which can be converted to an integer) with the
        minumum valid value for the parameter.

    maxi: str
        A string containing a function using n_col to convert number
        of columns to the maximum valid value for the paramter.

    Returns
    -------
    int_range_checker: function
        A function to input as a type to argparse to check if the value
        provided is in the right range
    '''
    # Based on solution from @georg-w stack overflow qu. 55324449
    mini = int(mini)
    try:
        maxi_val = eval(maxi)
    except TypeError:
        maxi_val = int(maxi)

    def int_range_checker(arg):
        try:
            f = int(arg)
        except ValueError:
            raise configargparse.ArgumentTypeError("Must be an integer")

        if f < mini or f > maxi_val:
            raise configargparse.ArgumentTypeError(
                "Must be in range [%s .. %s (%s)]" % (mini, maxi, maxi_val))
        return (f)

    return (int_range_checker)


def getParser():
    '''
    Builds a configargparse.ArgumentParser object with the CIAlign parameters

    Returns
    -------
    parser: configargparse.ArgumentParser
        ArgumentParser with the CIAlign parameters
    '''

    parser = configargparse.ArgumentParser(
             description='Clean and interpret a multiple sequence \
                          alignment', add_help=False)
    ci_dir = os.path.dirname(utilityFunctions.__file__)

    # Looks up the default values and minimum and maximum values for the
    # paramters associated with the cleaning functions in the text file
    # ranges.txt provided in the CIAlign code directory
    ranges = [line.strip().split("\t")
              for line in open("%s/ranges.txt" % ci_dir)]
    # Defaults
    defs = {x[0]: x[1] for x in ranges}
    # Minima
    minis = {x[0]: x[2] for x in ranges}
    # Maxima
    maxis = {x[0]: x[3] for x in ranges}

    # Seperate the required and optional paramters
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')

    # Files
    # not to confuse with inifile
    required.add("--infile", dest='infile', type=str,
                 help='Path to input alignment file in FASTA format')
    optional.add("--inifile", dest='inifile', type=str,
                 default=None,
                 help='Path to config file. Default: %(default)s',
                 is_config_file=True)
    optional.add("--outfile_stem", dest='outfile_stem', type=str,
                 default="CIAlign",
                 help="Prefix for output files, including the path to the \
                     output directory. Default: %(default)s")

    # Initial setup
    # Read the alignment temporarily just to find out how many columns there
    # are as for several of the cleaning functions the range of valid
    # parameters depends on this.

    tempargs = parser.parse_known_args()[0]
    if tempargs.infile:
        # Read the FASTA file into an array
        arr, nams = utilityFunctions.FastaToArray(tempargs.infile, None,
                                                  tempargs.outfile_stem)

        # Find the number of columns in the input alignment
        n_col = np.shape(arr)[1]
        # Remove the array from memory
        del arr
    else:
        # Gives a valid int value just for generating the --help text
        n_col = 100
    # parameter to run all functions without having to type them in
    optional.add("--all", dest="all_options",
                 action="store_true",
                 help="Use all available functions with default parameters.")

    # parameter to run all cleaning functions without having to type them in
    optional.add("--clean", dest="clean",
                 action="store_true",
                 help="Use all cleaning functions with default parameters.")

    # parameter to create all mini alignments without having to type them in
    optional.add("--visualise", dest="visualise",
                 action="store_true",
                 help="Plot all mini alignments with default parameters.")

    # parameter to run all interpreation functions except creating sequence
    # logos without having to type them in
    optional.add("--interpret", dest="interpret",
                 action="store_true",
                 help="Use all interpreting functions with default \
                 parameters.")

    # Runtime
    optional.add("--silent", dest='silent',
                 help="Do not print progress to the screen. \
                       Default: %(default)s",
                 action='store_true')

    # Crop Ends
    optional.add("--crop_ends", dest="crop_ends",
                 action="store_true",
                 help="Crop the ends of sequences if they are poorly aligned. \
                 Default: %(default)s")

    optional.add("--crop_ends_mingap_perc", dest='crop_ends_mingap_perc',
                 type=float_range(minis['crop_ends_mingap_perc'],
                                  maxis['crop_ends_mingap_perc']),
                 default=defs['crop_ends_mingap_perc'],
                 help="Minimum proportion of the sequence length (excluding \
                     gaps) that is the threshold for change in gap numbers. \
                     Default: %(default)s.",
                 metavar="(float, %s..%s)" % (minis['crop_ends_mingap_perc'],
                                              maxis['crop_ends_mingap_perc']))

    optional.add("--crop_ends_redefine_perc", dest='crop_ends_redefine_perc',
                 type=float_range(minis['crop_ends_redefine_perc'],
                                  maxis['crop_ends_redefine_perc']),
                 default=defs['crop_ends_redefine_perc'],
                 help="Proportion of the sequence length (excluding gaps) \
                       that is being checked for change in gap numbers to \
                       redefine start/end. Default: %(default)s",
                 metavar="(float, %s..%s)" % (
                     minis['crop_ends_redefine_perc'],
                     maxis['crop_ends_redefine_perc']))

    # Remove divergent sequences
    optional.add("--remove_divergent", dest="remove_divergent",
                 action="store_true",
                 help="Remove sequences with <= N proportion of positions at \
                       which the most common base / amino acid in the \
                       alignment is present. Default: %(default)s")

    optional.add("--remove_divergent_minperc", dest="remove_divergent_minperc",
                 default=defs['remove_divergent_minperc'],
                 type=float_range(minis['remove_divergent_minperc'],
                                  maxis['remove_divergent_minperc']),
                 help="Minimum proportion of positions which should be \
                       identical to the most common base / amino acid in \
                       order to be preserved. \
                       Default: %(default)s)",
                 metavar="(float, %s..%s)" % (
                     minis['remove_divergent_minperc'],
                     maxis['remove_divergent_minperc']))

    # # Remove Insertions
    optional.add("--remove_insertions", dest="remove_insertions",
                 action="store_true",
                 help="Remove insertions found in <= insertion_min_perc \
                 percent of sequences from the alignment. \
                 Default: %(default)s")

    optional.add("--insertion_min_size", dest="insertion_min_size",
                 type=int_range(minis['insertion_min_size'],
                                maxis['insertion_max_size'],
                                n_col),
                 default=defs['insertion_min_size'],
                 help="Only remove insertions >= this number of residues. \
                       Default: %(default)s.",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_min_size'],
                     maxis['insertion_min_size']))

    optional.add("--insertion_max_size", dest="insertion_max_size",
                 type=int_range(minis['insertion_max_size'],
                                maxis['insertion_max_size'],
                                n_col),
                 default=defs['insertion_max_size'],
                 help="Only remove insertions <= this number of residues. \
                       Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_max_size'],
                     maxis['insertion_max_size']))

    optional.add("--insertion_min_flank", dest="insertion_min_flank",
                 type=int_range(minis['insertion_min_flank'],
                                maxis['insertion_min_flank'],
                                n_col),
                 default=defs['insertion_min_flank'],
                 help="Minimum number of bases on either side of an insertion \
                       to classify it as an insertion.\
                       Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_min_flank'],
                     maxis['insertion_min_flank']))

    optional.add("--insertion_min_perc", dest="insertion_min_perc",
                 type=float_range(minis['insertion_min_perc'],
                                  maxis['insertion_min_perc']),
                 default=defs['insertion_min_perc'],
                 help="Remove insertions which are present in less than this \
                       proportion of sequences.\
                       Default: %(default)s",
                 metavar="(float, %s..%s)" % (
                     minis['insertion_min_perc'],
                     maxis['insertion_min_perc']))

    # Remove Short
    optional.add("--remove_short", dest="remove_short",
                 help="Remove sequences <= N bases / amino acids from the \
                       alignment. Default: %(default)s",
                 action="store_true")

    optional.add("--remove_min_length", dest="remove_min_length",
                 type=int_range(minis['remove_min_length'],
                                maxis['remove_min_length'],
                                n_col),
                 default=defs['remove_min_length'],
                 help="Sequences are removed if they are shorter than this \
                       minimum length, excluding gaps. Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['remove_min_length'],
                     maxis['remove_min_length']))

    # keep gap only
    optional.add("--keep_gaponly", dest="remove_gaponly",
                 action="store_false",
                 help="Keep gap only columns in the alignment. Default: \
                       %(default)s")

    # Consensus
    optional.add("--make_consensus", dest="make_consensus",
                 action="store_true",
                 help="Make a consensus sequence based on the cleaned \
                       alignment. Default: %(default)s")
    optional.add("--consensus_type", dest="consensus_type", type=str,
                 default="majority",
                 help="Type of consensus sequence to make - can be majority, \
                       to use the most common character at each position in \
                       the consensus, even if this is a gap, or \
                       majority_nongap, to use the most common non-gap \
                       character at each position. Default: %(default)s")
    optional.add("--consensus_keep_gaps", dest="consensus_keep_gaps",
                 action="store_true",
                 help="If there are gaps in the consensus (if majority_nongap \
                       is used as consensus_type), should these be included \
                       in the consensus (True) or should this position in \
                      the consensus be deleted (False). Default: %(default)s")
    optional.add("--consensus_name", dest="consensus_name",
                 type=str, default="consensus",
                 help="Name to use for the consensus sequence in the output \
                       fasta file. Default: %(default)s")

    # Mini Alignments
    optional.add("--plot_input", dest="plot_input",
                 action="store_true",
                 help="Plot a mini alignment - an image representing the \
                       input alignment. Default: %(default)s")
    optional.add("--plot_output", dest="plot_output",
                 action="store_true",
                 help="Plot a mini alignment, an image representing the \
                       output alignment. Default: %(default)s")
    optional.add("--plot_markup", dest="plot_markup",
                 action="store_true",
                 help="Draws the input alignment but with the columns and \
                       rows which have been removed by each function marked \
                       up in corresponding colours. Default: %(default)s")
    optional.add("--plot_dpi", dest="plot_dpi",
                 type=int, default=300,
                 help="DPI for mini alignments. Default: %(default)s")
    optional.add("--plot_format", dest="plot_format",
                 type=str, default='png',
                 help="Image format for mini alignments - can be png, svg, \
                       tiff or jpg. Default: %(default)s")
    optional.add("--plot_width", dest="plot_width",
                 type=int, default=5,
                 help="Mini alignment width in inches. Default: %(default)s")
    optional.add("--plot_height", dest="plot_height",
                 type=int, default=3,
                 help="Mini alignment height in inches. Default: %(default)s")
    optional.add("--plot_keep_numbers", dest="plot_keep_numbers",
                 action="store_true",
                 help="If specified, for mini alignments based on CIAlign \
                       output with <10 sequences (or if force_numbers \
                       is switched on) the rows will be labelled \
                       based on the input alignment, rather \
                       than renumbered")
    optional.add("--plot_force_numbers", dest="plot_force_numbers",
                 action="store_true",
                 help="Force all rows to be numbered on the mini alignments \
                 rather than labelling e.g. every 10th row for larger plots. \
                 Will cause labels to overlap on large plots")

    # Sequence logos
    optional.add("--make_sequence_logo", dest="make_sequence_logo",
                 action="store_true",
                 help="Draw a sequence logo. Default: %(default)s")
    optional.add("--sequence_logo_type", dest="sequence_logo_type",
                 type=str, default='bar',
                 help="Type of sequence logo - bar/text/both. \
                       Default: %(default)s")
    optional.add("--sequence_logo_dpi", dest="sequence_logo_dpi",
                 type=int, default=300,
                 help="DPI for sequence logo image. Default: %(default)s")
    optional.add("--sequence_logo_font", dest="sequence_logo_font",
                 type=str, default='monospace',
                 help="Font for text sequence logo. Default: %(default)s")
    optional.add("--sequence_logo_nt_per_row", dest='sequence_logo_nt_per_row',
                 type=int, default=50,
                 help="Number of bases / amino acids to show per row in the \
                       sequence logo, where the logo is too large to show on \
                       a single line. Default: %(default)s")
    optional.add("--sequence_logo_filetype", dest='sequence_logo_filetype',
                 type=str, default='png',
                 help="Image file type to use for the sequence logo - can be \
                       png, svg, tiff or jpg. Default: %(default)s")
    optional.add("--logo_start", dest="logo_start",
                 type=int, default=0,
                 help="Start position of sequence logo. Default: %(default)s")
    optional.add("--logo_end", dest="logo_end",
                 type=int, default=0,
                 help="End position of sequence logo. Default: %(default)s")
    optional.add("--list_fonts_only", dest='list_fonts_only',
                 action="store_true",
                 help="Make a swatch showing available fonts. \
                       Default: %(default)s")

    # Coverage
    optional.add("--plot_coverage_input", dest="plot_coverage_input",
                 action="store_true",
                 help="Plot the coverage of the input MSA. Default: \
                       %(default)s")
    optional.add("--plot_coverage_output", dest="plot_coverage_output",
                 action="store_true",
                 help="Plot the coverage of the output MSA. Default: \
                       %(default)s")
    optional.add("--plot_coverage_dpi", dest="plot_coverage_dpi",
                 type=int, default=300,
                 help="DPI for coverage plot. Default: %(default)s")
    optional.add("--plot_coverage_height", dest="plot_coverage_height",
                 type=int, default=3,
                 help="Height for coverage plot (inches). Default: \
                       %(default)s")
    optional.add("--plot_coverage_width", dest="plot_coverage_width",
                 type=int, default=5,
                 help="Width for coverage plot (inches). Default: \
                       %(default)s")
    optional.add("--plot_coverage_colour", dest="plot_coverage_colour",
                 type=str, default='#007bf5',
                 help="Colour for coverage plot (hex code or name). \
                       Default: %(default)s")
    optional.add("--plot_coverage_filetype", dest="plot_coverage_filetype",
                 type=str, default='png',
                 help="File type for coverage plot (png, svg, tiff, jpg). \
                       Default: %(default)s")

    # Similarity Matrix
    optional.add("--make_similarity_matrix_input", dest="make_simmatrix_input",
                 action="store_true",
                 help="Make a similarity matrix for the input alignment. \
                       Default: %(default)s")
    optional.add("--make_similarity_matrix_output",
                 dest="make_simmatrix_output",
                 action="store_true",
                 help="Make a similarity matrix for the output alignment. \
                       Default: %(default)s")
    optional.add("--make_simmatrix_dp", dest="make_simmatrix_dp",
                 type=int, default=4,
                 help="Number of decimal places to display in the similarity \
                       matrix output file. Default: %(default)s")
    optional.add("--make_simmatrix_minoverlap",
                 dest="make_simmatrix_minoverlap",
                 type=int, default=1,
                 help="Minimum overlap between two sequences to have non-zero \
                       similarity in the similarity matrix. \
                       Default: %(default)s")
    optional.add("--make_simmatrix_keepgaps", dest="make_simmatrix_keepgaps",
                 type=int, default=0,
                 help="Include positions with gaps in either or both \
                       sequences in the similarity matrix calculation. \
                       Can be 0 - exclude positions which are gaps in either \
                       or both sequences, 1 - exclude positions which are \
                       gaps in both sequences, 2 - consider all positions \
                       regardless of gaps. Default: %(default)s")

    # Unalign function
    optional.add("--unalign_input", dest="unalign_input",
                 action="store_true", default=False,
                 help="Generate a copy of the input alignment with no gaps. \
                       Default: %(default)s")
    optional.add("--unalign_output", dest="unalign_output",
                 action="store_true", default=False,
                 help="Generate a copy of the cleaned alignment with no \
                     gaps. Default: %(default)s")

    # Replace Us by Ts function
    optional.add("--replace_input", dest="replace_input", action="store_true",
                 default=False,
                 help="Replaces all Us by Ts in input alignment. \
                     Default: %(default)s")
    optional.add("--replace_output", dest="replace_output",
                 action="store_true", default=False,
                 help="Replaces all Us by Ts in output alignment. \
                     Default: %(default)s")

    # Help function
    optional.add('-h', '--help', action='help',
                 default=configargparse.SUPPRESS,
                 help='Show all available parameters with an explanation.')

    # Version function
    optional.add('-v', '--version', action='version',
                 version=__version__,
                 default=configargparse.SUPPRESS,
                 help='Show the current version.')
    return (parser)
