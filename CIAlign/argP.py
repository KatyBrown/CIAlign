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


def float_range(mini, maxi, default):
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
    default = float(default)

    def float_range_checker(arg):
        try:
            f = float(arg)
        except ValueError:
            raise configargparse.ArgumentTypeError(
                "Must be a floating point number")
        if (f < mini or f > maxi) and not f == default:
            raise configargparse.ArgumentTypeError(
                "Must be in range [%s .. %s]" % (mini, maxi))
        return (f)

    return (float_range_checker)


def int_range(mini, maxi, n_col, default):
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
    try:
        default_val = eval(default)
    except TypeError:
        default_val = int(default)

    def int_range_checker(arg):
        try:
            f = int(arg)
        except ValueError:
            raise configargparse.ArgumentTypeError("Must be an integer")

        if (f < mini or f > maxi_val) and not f == default_val:
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
                 help='Path to input alignment file in FASTA format. \
                       Required')
    optional.add("--inifile", dest='inifile', type=str,
                 default=None,  metavar="(string)",
                 help='Path to config file. Default: %(default)s',
                 is_config_file=True)
    optional.add("--outfile_stem", dest='outfile_stem', type=str,
                 default="CIAlign", metavar="(string)",
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
                 help="""Use all available functions, with default parameters \
                         unless others are specified. \
                         Default: %(default)s""")

    # parameter to run all cleaning functions without having to type them in
    optional.add("--clean", dest="clean",
                 action="store_true",
                 help="""Use all cleaning functions, with default parameters \
                         unless others are specified. \
                         Default: %(default)s""")

    # parameter to create all mini alignments without having to type them in
    optional.add("--visualise", dest="visualise",
                 action="store_true",
                 help="""Plot all mini alignments, with default parameters \
                         unless others are specified. \
                         Default: %(default)s""")

    # parameter to run all interpreation functions except creating sequence
    # logos without having to type them in
    optional.add("--interpret", dest="interpret",
                 action="store_true",
                 help="""Use all interpretation functions, with default
                         parameters unless others are specified. \
                         Default: %(default)s""")

    # Runtime
    optional.add("--silent", dest='silent',
                 help="Do not print progress to the screen. \
                       Default: %(default)s",
                 action='store_true')

    # Remove divergent sequences
    optional.add("--remove_divergent", dest="remove_divergent",
                 action="store_true",
                 help="Remove sequences with <= N proportion of positions at \
                       which the most common base / amino acid in the \
                       alignment is present. Default: %(default)s")

    optional.add("--remove_divergent_minperc", dest="remove_divergent_minperc",
                 default=defs['remove_divergent_minperc'],
                 type=float_range(minis['remove_divergent_minperc'],
                                  maxis['remove_divergent_minperc'],
                                  defs['remove_divergent_minperc']),
                 help="Minimum proportion of positions which should be \
                       identical to the most common base / amino acid in \
                       order to be preserved. \
                       Default: %(default)s",
                 metavar="(float, %s..%s)" % (
                     minis['remove_divergent_minperc'],
                     maxis['remove_divergent_minperc']))

    optional.add("--remove_divergent_retain", dest="retain_seqs_rd",
                 action="append", default=None, type=str,
                 metavar="(string)",
                 help="""Do not remove the sequence with this name when \
                         running the remove_divergent function. \
                         Sequence names must exactly match the FASTA infile. \
                         Can be specified \
                         multiple times. Default: %(default)s""")

    optional.add("--remove_divergent_retain_str", dest="retain_seqs_rdS",
                 action="append", default=None, type=str,
                 metavar="(string)",
                 help="""Do not remove the sequences with names containing \
                         this word (character string) when \
                         running the remove_divergent function. \
                         Case sensitive. \
                         Default: %(default)s""")

    optional.add("--remove_divergent_retain_list", dest="retain_seqs_rdL",
                 type=str, default=None,
                 metavar="(string)",
                 help="""Do not remove the sequences listed in this file when \
                         running the remove_divergent function. \
                         Sequence names must exactly match the FASTA infile. \
                         Default: %(default)s""")

    # # Remove Insertions
    optional.add("--remove_insertions", dest="remove_insertions",
                 action="store_true",
                 help="Remove insertions found in <= insertion_min_perc \
                 percent of sequences from the alignment. \
                 Default: %(default)s")

    optional.add("--insertion_min_size", dest="insertion_min_size",
                 type=int_range(minis['insertion_min_size'],
                                maxis['insertion_min_size'],
                                n_col, defs['insertion_min_size']),
                 default=defs['insertion_min_size'],
                 help="Only remove insertions >= this number of residues. \
                       Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_min_size'],
                     maxis['insertion_min_size']))

    optional.add("--insertion_max_size", dest="insertion_max_size",
                 type=int_range(minis['insertion_max_size'],
                                maxis['insertion_max_size'],
                                n_col,
                                defs['insertion_max_size']),
                 default=defs['insertion_max_size'],
                 help="Only remove insertions <= this number of residues. \
                       Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_max_size'],
                     maxis['insertion_max_size']))

    optional.add("--insertion_min_flank", dest="insertion_min_flank",
                 type=int_range(minis['insertion_min_flank'],
                                maxis['insertion_min_flank'],
                                n_col,
                                defs['insertion_min_flank']),
                 default=defs['insertion_min_flank'],
                 help="Minimum number of bases on either side of an insertion \
                       to classify it as an insertion.\
                       Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['insertion_min_flank'],
                     maxis['insertion_min_flank']))

    optional.add("--insertion_min_perc", dest="insertion_min_perc",
                 type=float_range(minis['insertion_min_perc'],
                                  maxis['insertion_min_perc'],
                                  defs['insertion_min_perc']),
                 default=defs['insertion_min_perc'],
                 help="Remove insertions which are present in less than this \
                       proportion of sequences.\
                       Default: %(default)s",
                 metavar="(float, %s..%s)" % (
                     minis['insertion_min_perc'],
                     maxis['insertion_min_perc']))

    # Crop Ends
    optional.add("--crop_ends", dest="crop_ends",
                 action="store_true",
                 help="Crop the ends of sequences if they are poorly aligned. \
                 Default: %(default)s")

    optional.add("--crop_ends_mingap_perc", dest='crop_ends_mingap_perc',
                 type=float_range(minis['crop_ends_mingap_perc'],
                                  maxis['crop_ends_mingap_perc'],
                                  defs['crop_ends_mingap_perc']),
                 default=defs['crop_ends_mingap_perc'],
                 help="Minimum proportion of the sequence length (excluding \
                     gaps) that is the threshold for change in gap numbers. \
                     Default: %(default)s",
                 metavar="(float, %s..%s)" % (minis['crop_ends_mingap_perc'],
                                              maxis['crop_ends_mingap_perc']))

    optional.add("--crop_ends_redefine_perc", dest='crop_ends_redefine_perc',
                 type=float_range(minis['crop_ends_redefine_perc'],
                                  maxis['crop_ends_redefine_perc'],
                                  defs['crop_ends_redefine_perc']),
                 default=defs['crop_ends_redefine_perc'],
                 help="Proportion of the sequence length (excluding gaps) \
                       that is being checked for change in gap numbers to \
                       redefine start/end. Default: %(default)s",
                 metavar="(float, %s..%s)" % (
                     minis['crop_ends_redefine_perc'],
                     maxis['crop_ends_redefine_perc']))

    optional.add("--crop_ends_retain", dest="retain_seqs_ce",
                 action="append", default=None, metavar="(string)",
                 help="""Do not crop the sequence with this name when \
                         running the crop_ends function. Can be specified \
                         multiple times. Default: %(default)s""")

    optional.add("--crop_ends_retain_str", dest="retain_seqs_ceS",
                 action="append", default=None, type=str,
                 metavar="(string)",
                 help="""Do not crop sequences with names containing \
                         this word (character string) when \
                         running the crop_ends function. \
                         Case sensitive. \
                         Default: %(default)s""")

    optional.add("--crop_ends_retain_list", dest="retain_seqs_ceL",
                 type=str, default=None,
                 metavar="(string)",
                 help="""Do not crop the sequences listed in this file when \
                         running the crop_ends function. \
                         Sequence names must exactly match the FASTA infile. \
                         Default: %(default)s""")

    # Remove Short
    optional.add("--remove_short", dest="remove_short",
                 help="Remove sequences <= N bases / amino acids from the \
                       alignment. Default: %(default)s",
                 action="store_true")

    optional.add("--remove_min_length", dest="remove_min_length",
                 type=int_range(minis['remove_min_length'],
                                maxis['remove_min_length'],
                                n_col,
                                defs['remove_min_length']),
                 default=defs['remove_min_length'],
                 help="Sequences are removed if they are shorter than this \
                       minimum length, excluding gaps. Default: %(default)s",
                 metavar="(int, %s..%s)" % (
                     minis['remove_min_length'],
                     maxis['remove_min_length']))
    optional.add("--remove_short_retain", dest="retain_seqs_rs",
                 action="append", default=None, metavar="(string)",
                 help="""Do not remove the sequence with this name when \
                         running the remove_divergent function.
                         Sequence names must exactly match the FASTA infile. \
                         Can be specified multiple times. \
                         Default: %(default)s""")

    optional.add("--remove_short_retain_str", dest="retain_seqs_rsS",
                 action="append", default=None, type=str, metavar="(string)",
                 help="""Do not remove the sequences with names containing \
                         this word (character string) when \
                         running the remove_short function. \
                         Case sensitive. \
                         Default: %(default)s""")

    optional.add("--remove_short_retain_list", dest="retain_seqs_rsL",
                 type=str, default=None, metavar="(string)",
                 help="""Do not remove the sequences listed in this file when \
                         running the remove_short function. \
                         Sequence names must exactly match the FASTA infile. \
                         Default: %(default)s""")

    # keep gap only
    optional.add("--keep_gaponly", dest="remove_gaponly",
                 action="store_false",
                 help="Keep gap only columns in the alignment. Default: \
                       %(default)s")

    # Crop Divergent
    optional.add("--crop_divergent", dest="crop_divergent",
                 help="""Crop ends of sequences which are highly
                         divergent. Default: %(default)s""",
                 action="store_true")

    optional.add("--crop_divergent_min_prop_ident",
                 dest="divergent_min_prop_ident",
                 type=float_range(minis['divergent_min_prop_ident'],
                                  maxis['divergent_min_prop_ident'],
                                  defs['divergent_min_prop_ident']),
                 default=defs['divergent_min_prop_ident'],
                 help="""The minimum proportion of sequences which should \
                         have the same residue in each column for crop \
                         divergent. Default: %(default)s""",
                 metavar="(float, %s..%s)" % (
                     minis['divergent_min_prop_ident'],
                     maxis['divergent_min_prop_ident']))

    optional.add("--crop_divergent_min_prop_nongap",
                 dest="divergent_min_prop_nongap",

                 type=float_range(minis['divergent_min_prop_nongap'],
                                  maxis['divergent_min_prop_nongap'],
                                  defs['divergent_min_prop_nongap']),
                 default=defs['divergent_min_prop_nongap'],
                 help="""The minimum proportion of sequences which should \
                         have the non-gap residues in each column \
                         for crop divergent. Default: %(default)s""",
                 metavar="(float, %s..%s)" % (
                     minis['divergent_min_prop_nongap'],
                     maxis['divergent_min_prop_nongap']))

    optional.add("--crop_divergent_buffer_size",
                 dest="divergent_buffer_size",
                 type=int_range(minis['divergent_buffer_size'],
                                maxis['divergent_buffer_size'],
                                n_col,
                                defs['divergent_buffer_size']),
                 default=defs['divergent_buffer_size'],
                 help="""The number of consecutive columns which should meet \
                         the min_prop_ident and min_prop_nongap criteria \
                         to pass filtering by crop_divergent. \
                         Default: %(default)s""",
                 metavar="(int, %s..%s)" % (
                     minis['divergent_buffer_size'],
                     maxis['divergent_buffer_size']))

    # Retain
    optional.add("--retain", dest="retain_seqs",
                 action="append", default=None, type=str, metavar="(string)",
                 help="""Do not remove the sequence with this name when \
                         running any rowwise function \
                         (currently remove_divergent and crop_ends). \
                         Sequence names must exactly match the FASTA infile. \
                         Can be specified multiple times. \
                         Default: %(default)s""")

    optional.add("--retain_str", dest="retain_seqsS",
                 action="append", default=None, metavar="(string)",
                 help="""Do not remove the sequences with names containing \
                         this word (character string) when \
                         running any rowwise function \
                         (currently remove_divergent and crop_ends). \
                         Case sensitive. \
                         Default: %(default)s""")

    optional.add("--retain_list", dest="retain_seqsL",
                 type=str, default=None, metavar="(string)",
                 help="""Do not remove the sequences listed in this file when \
                         running any rowwise function \
                         (currently remove_divergent, remove_short and\
                          crop_ends). \
                         Sequence names must exactly match the FASTA infile. \
                         Default: %(default)s""")

    # Consensus
    optional.add("--make_consensus", dest="make_consensus",
                 action="store_true",
                 help="Make a consensus sequence based on the cleaned \
                       alignment. Default: %(default)s")

    optional.add("--consensus_type", dest="consensus_type", type=str,
                 default="majority", metavar="(string)",
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
                 type=str, default="consensus", metavar="(string)",
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
                 type=int, default=300, metavar="(int)",
                 help="DPI for mini alignments. Default: %(default)s")

    optional.add("--plot_format", dest="plot_format",
                 type=str, default='png', metavar="(string)",
                 help="Image format for mini alignments - can be png, svg, \
                       tiff or jpg. Default: %(default)s")

    optional.add("--plot_width", dest="plot_width",
                 type=int, default=5, metavar="(int)",
                 help="Mini alignment width in inches. Default: %(default)s")

    optional.add("--plot_height", dest="plot_height",
                 type=int, default=3, metavar="(int)",
                 help="Mini alignment height in inches. Default: %(default)s")

    optional.add("--plot_keep_numbers", dest="plot_keep_numbers",
                 action="store_true",
                 help="If specified, for mini alignments based on CIAlign \
                       output with <10 sequences (or if force_numbers \
                       is switched on) the rows will be labelled \
                       based on the input alignment, rather \
                       than renumbered. Default: %(default)s")

    optional.add("--plot_force_numbers", dest="plot_force_numbers",
                 action="store_true",
                 help="Force all rows to be numbered on the mini alignments \
                 rather than labelling e.g. every 10th row for larger plots. \
                 Will cause labels to overlap on large plots. \
                 Default: %(default)s")

    # Colours
    optional.add("--palette", dest="palette", type=str,
                 default="CBS", metavar="(str",
                 help="Colour palette. Currently implemented \
                       CBS (colour blind safe) or bright. \
                       Default: %(default)s")

    # Sequence logos
    optional.add("--make_sequence_logo", dest="make_sequence_logo",
                 action="store_true",
                 help="Draw a sequence logo. Default: %(default)s")

    optional.add("--sequence_logo_type", dest="sequence_logo_type",
                 type=str, default='bar', metavar="(string)",
                 help="Type of sequence logo - bar/text/both. \
                       Default: %(default)s")

    optional.add("--sequence_logo_dpi", dest="sequence_logo_dpi",
                 type=int, default=300, metavar="(int)",
                 help="DPI for sequence logo image. Default: %(default)s")

    optional.add("--sequence_logo_font", dest="sequence_logo_font",
                 type=str, default='monospace', metavar="(string)",
                 help="Font for text sequence logo. Default: %(default)s")

    optional.add("--sequence_logo_nt_per_row", dest='sequence_logo_nt_per_row',
                 type=int, default=50, metavar="(int)",
                 help="Number of bases / amino acids to show per row in the \
                       sequence logo, where the logo is too large to show on \
                       a single line. Default: %(default)s")

    optional.add("--sequence_logo_filetype", dest='sequence_logo_filetype',
                 type=str, default='png', metavar="(string)",
                 help="Image file type to use for the sequence logo - can be \
                       png, svg, tiff or jpg. Default: %(default)s")

    optional.add("--logo_start", dest="logo_start",
                 type=int, default=0, metavar="(int)",
                 help="Start position of sequence logo. Default: %(default)s")

    optional.add("--logo_end", dest="logo_end",
                 type=int, default=0, metavar="(int",
                 help="End position of sequence logo. Default: %(default)s")

    optional.add("--list_fonts_only", dest='list_fonts_only',
                 action="store_true",
                 help="Make a swatch showing available fonts. \
                       Default: %(default)s")

    # PWM function
    optional.add("--pwm_input", dest="pwm_input",
                 action="store_true", default=False,
                 help="Generate a position frequency matrix, position \
                       probability matrix and position weight matrix based \
                       on the input alignment. Default: %(default)s")

    optional.add("--pwm_output", dest="pwm_output",
                 action="store_true", default=False,
                 help="Generate a position frequency matrix, position \
                       probability matrix and position weight matrix based \
                       on the output alignment. Default: %(default)s")

    optional.add("--pwm_start", dest="pwm_start",
                 type=int, default=None, metavar="(int)",
                 help="Start column of the PWM. Default: %(default)s")

    optional.add("--pwm_end", dest="pwm_end",
                 type=int, default=None, metavar="(int",
                 help="End column of the PWM. Default: %(default)s")

    optional.add("--pwm_freqtype", dest="pwm_freqtype",
                 type=str, default="equal", metavar="(str",
                 help="Type of background frequency matrix to use when \
                       generating the PWM. Should be 'equal', 'calc', 'calc2' \
                       or user. 'equal', assume all residues are equally \
                       common, 'calc', frequency is calculated using the PFM, \
                       'calc2', frequency is calculated using the full \
                       alignment (same as calc if pwm_start and pwm_end are \
                       not specified). Default: %(default)s")

    optional.add("--pwm_alphatype", dest="pwm_alphatype",
                 type=str, default="calc", metavar="(str",
                 help="Alpha value to use as a pseudocount to avoid zero \
                       values in the PPM. Should be 'calc' or 'user'. \
                       If alphatype is 'calc', alpha is calculated as \
                       frequency(base) * (square root(n rows in alignment)), \
                       as described in Dave Tang's blog here: \
                       https://davetang.org/muse/2013/10/01/\
                       position-weight-matrix/, \
                       which recreates the method used in \
                       doi.org/10.1038/nrg1315. If alpha type is 'user' \
                       the user provides the value of alpha as pwm_alphatype. \
                       To run without pseudocounts set pwm_alphatype as user \
                       and pwm_alphaval as 0. Default: %(default)s")

    optional.add("--pwm_alphaval", dest="pwm_alphaval",
                 type=float, default=1.0, metavar="(int",
                 help="User defined value of the alpha parameter to use as a \
                       pseudocount in the PPM. Default: %(default)s")

    optional.add("--pwm_output_blamm", dest="pwm_output_blamm",
                 action="store_true", default=False,
                 help="Output PPM formatted for BLAMM software \
                       https://github.com/biointec/blamm. \
                       Default: %(default)s")

    optional.add("--pwm_output_meme", dest="pwm_output_meme",
                 action="store_true", default=False,
                 help="Output PPM formatted for MEME software \
                       https://meme-suite.org/meme \
                       Default: %(default)s")

    # Plots
    optional.add("--plot_stats_input", dest="plot_stats_input",
                 action="store_true",
                 help="Plot statistics about the input MSA. Default: \
                       %(default)s")

    optional.add("--plot_stats_output", dest="plot_stats_output",
                 action="store_true",
                 help="Plot statistics about the output MSA. Default: \
                       %(default)s")

    optional.add("--plot_stats_dpi", dest="plot_stats_dpi",
                 type=int, default=300, metavar="(int)",
                 help="DPI for coverage plot. Default: %(default)s")

    optional.add("--plot_stats_height", dest="plot_stats_height",
                 type=int, default=3, metavar="(int)",
                 help="Height for statistics plots (inches). Default: \
                       %(default)s")

    optional.add("--plot_stats_width", dest="plot_stats_width",
                 type=int, default=5, metavar="(int)",
                 help="Width for statistics plots (inches). Default: \
                       %(default)s")

    optional.add("--plot_stats_colour", dest="plot_stats_colour",
                 type=str, default='#007bf5', metavar="(string)",
                 help="Colour for statistics plots (hex code or name). \
                       Default: %(default)s")

    optional.add("--plot_stats_filetype", dest="plot_stats_filetype",
                 type=str, default='png', metavar="(string)",
                 help="File type for statistics plots (png, svg, tiff, jpg). \
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
                 type=int, default=4, metavar="(int)",
                 help="Number of decimal places to display in the similarity \
                       matrix output file. Default: %(default)s")

    optional.add("--make_simmatrix_minoverlap",
                 dest="make_simmatrix_minoverlap",
                 type=int, default=1, metavar="(int)",
                 help="Minimum overlap between two sequences to have non-zero \
                       similarity in the similarity matrix. \
                       Default: %(default)s")

    optional.add("--make_simmatrix_keepgaps", dest="make_simmatrix_keepgaps",
                 type=int, default=0, metavar="(int)",
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
    optional.add("--replace_input_ut", dest="replace_input_ut",
                 action="store_true",
                 default=False,
                 help="Replaces all Us by Ts in input alignment. \
                     Default: %(default)s")

    optional.add("--replace_output_ut", dest="replace_output_ut",
                 action="store_true", default=False,
                 help="Replaces all Us by Ts in output alignment. \
                     Default: %(default)s")
    # Replace Ts by Us function
    optional.add("--replace_input_tu", dest="replace_input_tu",
                 action="store_true",
                 default=False,
                 help="Replaces all Ts by Us in input alignment. \
                     Default: %(default)s")

    optional.add("--replace_output_tu", dest="replace_output_tu",
                 action="store_true", default=False,
                 help="Replaces all Ts by Us in output alignment. \
                     Default: %(default)s")

    # Section
    optional.add("--get_section", dest="get_section",
                 action="store_true", default=False,
                 help="Retrieve and process a section of the alignment, \
                       requires the \
                       section_start and section_end parameters. All \
                       logging is relative to the original start position. \
                       Default: %(default)s")

    optional.add("--section_start", dest="section_start",
                 type=int, default=None, metavar="(int",
                 help="Start position (column) for a section of the alignment \
                 to be isolated. 0-based - the first column is column 0. \
                 Default: %(default)s")

    optional.add("--section_end", dest="section_end",
                 type=int, default=None, metavar="(int",
                 help="End position (column) for a section of the alignment \
                 to be isolated. 0-based - the first column is column 0. \
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
