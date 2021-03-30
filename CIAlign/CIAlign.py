#! /usr/bin/env python

import logging
import configargparse
import os.path
import numpy as np
import copy

try:
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.parsingFunctions as parsingFunctions
    import CIAlign.miniAlignments as miniAlignments
    import CIAlign.similarityMatrix as similarityMatrix
    import CIAlign.consensusSeq as consensusSeq
    from CIAlign._version import __version__
except ImportError:
    import utilityFunctions
    import parsingFunctions
    import miniAlignments
    import similarityMatrix
    import consensusSeq
    from _version import __version__


def main():

    parser = configargparse.ArgumentParser(
                description='Clean and interpret a multiple sequence \
                             alignment',
                add_help=False)
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

    # parameter to run all functions without having to type them in
    optional.add("--all", dest="all_options",
                 action="store_true",
                 help="Use all available functions with default parameters.")

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
                 type=float, default=0.05,
                 help="Minimum proportion of the sequence length (excluding \
                     gaps) that is the threshold for change in gap numbers. \
                     Default: %(default)s")
    optional.add("--crop_ends_redefine_perc", dest='crop_ends_redefine_perc',
                 type=float, default=0.1,
                 help="Proportion of the sequence length (excluding gaps) \
                       that is being checked for change in gap numbers to \
                       redefine start/end. Default: %(default)s")

    # Remove divergent sequences
    optional.add("--remove_divergent", dest="remove_divergent",
                 action="store_true",
                 help="Remove sequences with <= N proportion of positions at \
                       which the most common base / amino acid in the \
                       alignment is present. Default: %(default)s")
    optional.add("--remove_divergent_minperc", dest="remove_divergent_minperc",
                 type=float, default=0.65,
                 help="Minimum proportion of positions which should be \
                       identical to the most common base / amino acid in \
                       order to be preserved. Default: %(default)s")

    # # Remove Insertions
    optional.add("--remove_insertions", dest="remove_insertions",
                 action="store_true",
                 help="Remove insertions found in <= 50 percent of sequences \
                       from the alignment. Default: %(default)s")
    optional.add("--insertion_min_size", dest="insertion_min_size",
                 type=int, default=3,
                 help="Only remove insertions >= this number of residues. \
                       Default: %(default)s")
    optional.add("--insertion_max_size", dest="insertion_max_size",
                 type=int, default=200,
                 help="Only remove insertions <= this number of residues. \
                       Default: %(default)s")
    optional.add("--insertion_min_flank", dest="insertion_min_flank",
                 type=int, default=5,
                 help="Minimum number of bases on either side of an insertion \
                       to classify it as an insertion. Default: %(default)s")

    # Remove Short
    optional.add("--remove_short", dest="remove_short",
                 help="Remove sequences <= N bases / amino acids from the \
                       alignment. Default: %(default)s",
                 action="store_true")
    optional.add("--remove_min_length", dest="remove_min_length",
                 type=int, default=50,
                 help="Sequences are removed if they are shorter than this \
                       minimum length, excluding gaps. Default: %(default)s")

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

    args = parser.parse_args()

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s_log.txt" % args.outfile_stem
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)

    rmfile = "%s_removed.txt" % args.outfile_stem
    outfile = "%s_cleaned.fasta" % (args.outfile_stem)

    # add the handlers to the logger
    log.addHandler(handler)

    log.info("\nInitial parameters:\n%s" % str(parser.format_values()))

    if args.list_fonts_only:
        print("Listing fonts")
        out = "%s_fonts.png" % args.outfile_stem
        utilityFunctions.listFonts(out)
        exit()

    # convert the input fasta file into an array and make a list of
    # sequence names so the order can be maintained
    # XXXXX comes from the ini file
    if not args.infile or args.infile == 'XXXXX':
        print("Input alignment must be provided.")
        exit()

    # check if input file exists and is a file
    if not os.path.isfile(args.infile):
        print("Error! Your input alignmnent path could not be found.")
        exit()

    arr, nams = utilityFunctions.FastaToArray(args.infile, log, args.outfile_stem)

    # check if at least names are unique
    if len(nams) > len(set(nams)):
        print("Error! Your input alignmnent has duplicate names!")
        exit()

    # output file to store memory and time
    cmt_file = args.outfile_stem + "_cmt"

    cleaningArgs = [args.remove_insertions,
                    args.crop_ends,
                    args.remove_divergent]
    # check numbers of sequences first

    if len(arr) < 3 and any(cleaningArgs):
        # when less than three sequences, stop
        print("You need at least three sequences in your MSA to run \
               remove_insertions, crop_ends or remove_divergent")
        exit()
    elif len(arr) < 2:
        print("You need at least two sequences in your MSA")
        exit()

    # store a copy of the original array
    orig_arr = copy.copy(arr)
    orig_nams = copy.copy(nams)

    # make a dictionary to store the changes made
    markupdict = dict()

    # remember positions relative to original alignment
    relativePositions = list(range(0, len(orig_arr[0])))

    removed_seqs = set()
    removed_cols = set()
    removed_positions = dict()

    # detect if the sequence is amino acids or nucleotides
    typ = utilityFunctions.seqType(arr)

    if typ == 'aa':
        log.info("Amino acid alignment detected")
    else:
        log.info("Nucleotide alignment detected")

    rmfile = "%s_removed.txt" % args.outfile_stem
    outfile = "%s_cleaned.fasta" % (args.outfile_stem)
    reset_rmfile = open(rmfile, "w")
    reset_rmfile.close()
    utilityFunctions.checkArrLength(arr, log)

    if args.remove_divergent or args.all_options:
        log.info("Removing divergent sequences")
        if not args.silent:
            print("Removing divergent sequences")
        minperc = args.remove_divergent_minperc
        arr, r = parsingFunctions.removeDivergent(arr, nams,
                                                  rmfile, log,
                                                  minperc)

        markupdict['remove_divergent'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)
        utilityFunctions.checkArrLength(arr, log)

    if (args.remove_divergent and args.remove_gaponly) or args.all_options:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")
        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A

        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_insertions or args.all_options:
        log.info("Removing insertions")
        if not args.silent:
            print("Removing insertions")

        A = parsingFunctions.removeInsertions(arr,
                                              relativePositions,
                                              rmfile,
                                              log,
                                              args.insertion_min_size,
                                              args.insertion_max_size,
                                              args.insertion_min_flank)

        arr, r, relativePositions = A
        markupdict['remove_insertions'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if (args.remove_insertions and args.remove_gaponly) or args.all_options:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.crop_ends or args.all_options:
        # doesn't remove any whole columns or rows
        log.info("Cropping ends")
        if not args.silent:
            print("Cropping ends")
        arr, r = parsingFunctions.cropEnds(arr, nams, relativePositions,
                                           rmfile,
                                           log, args.crop_ends_mingap_perc,
                                           args.crop_ends_redefine_perc)

        markupdict['crop_ends'] = r
        removed_positions.update(r)
        utilityFunctions.checkArrLength(arr, log)

    if (args.crop_ends and args.remove_gaponly) or args.all_options:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_short or args.all_options:
        log.info("Removing short sequences")
        if not args.silent:
            print("Removing short sequences")
        arr, r = parsingFunctions.removeTooShort(arr, nams, rmfile, log,
                                                 args.remove_min_length)

        markupdict['remove_short'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)
        utilityFunctions.checkArrLength(arr, log)

    if (args.remove_short and args.remove_gaponly) or args.all_options:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_gaponly and not (args.all_options or
                                    args.remove_divergent or
                                    args.remove_insertions or
                                    args.crop_ends or
                                    args.remove_short):
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.make_simmatrix_input or args.all_options:
        log.info("Building similarity matrix for input alignment")
        if not args.silent:
            print("Building similarity matrix for input alignment")
        outf = "%s_input_similarity.tsv" % (args.outfile_stem)
        minoverlap = args.make_simmatrix_minoverlap
        keepgaps = args.make_simmatrix_keepgaps
        dp = args.make_simmatrix_dp
        similarityMatrix.calculateSimilarityMatrix(orig_arr,
                                                   orig_nams,
                                                   minoverlap=minoverlap,
                                                   keepgaps=keepgaps,
                                                   outfile=outf, dp=dp)

    if args.make_simmatrix_output or args.all_options:
        log.info("Building similarity matrix for output alignment")
        if not args.silent:
            print("Building similarity matrix for output alignment")
        outf = "%s_output_similarity.tsv" % (args.outfile_stem)
        minoverlap = args.make_simmatrix_minoverlap
        keepgaps = args.make_simmatrix_keepgaps
        dp = args.make_simmatrix_dp
        similarityMatrix.calculateSimilarityMatrix(arr,
                                                   nams,
                                                   minoverlap=minoverlap,
                                                   keepgaps=keepgaps,
                                                   outfile=outf, dp=dp)

    if args.plot_input or args.all_options:
        log.info("Plotting mini alignment for input")
        if not args.silent:
            print("Plotting mini alignment for input")
        outf = "%s_input.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, log, cmt_file,
                                         outf, typ, args.plot_dpi,
                                         False, args.plot_width,
                                         args.plot_height)

    if args.plot_output or args.all_options:
        log.info("Plotting mini alignment for output")
        if not args.silent:
            print("Plotting mini alignment for output")
        outf = "%s_output.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(arr, nams, log,
                                         outf, typ,
                                         args.plot_dpi,
                                         False,
                                         args.plot_width,
                                         args.plot_height)

    if args.plot_markup or args.all_options:
        log.info("Plotting mini alignment with markup")
        if not args.silent:
            print("Plotting mini alignment with markup")
        outf = "%s_markup.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, orig_nams, log,
                                         outf, typ,
                                         args.plot_dpi,
                                         False,
                                         args.plot_width, args.plot_height,
                                         markup=True, markupdict=markupdict)

    if args.make_consensus or args.all_options:
        log.info("Building consensus sequence")
        if not args.silent:
            print("Building consensus sequence")
        cons, coverage = consensusSeq.findConsensus(arr,
                                                    log, args.consensus_type)
        consarr = np.array(cons)
        arr_plus_cons = np.row_stack((arr, consarr))
        cons = "".join(cons)
        # Remove the gaps from the consensus if this option is specified
        if not args.consensus_keep_gaps:
            cons = cons.replace("-", "")
        out = open("%s_consensus.fasta" % args.outfile_stem, "w")
        out.write(">%s\n%s\n" % (args.consensus_name, cons))
        out.close()
        outf = "%s_with_consensus.fasta" % args.outfile_stem
        utilityFunctions.writeOutfile(outf, arr_plus_cons,
                                      nams + [args.consensus_name],
                                      removed_seqs)

    if args.plot_coverage_input or args.all_options:
        log.info("Plotting coverage for input")
        if not args.silent:
            print("Plotting coverage for input")
        coverage_file = "%s_input_coverage.%s" % (args.outfile_stem,
                                                  args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(orig_arr,
                                                     args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)

    if args.plot_coverage_output or args.all_options:
        if not args.silent:
            print("Plotting coverage for output")
        log.info("Plotting coverage for output")
        coverage_file = "%s_output_coverage.%s" % (args.outfile_stem,
                                                   args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(arr, args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)

    if args.make_sequence_logo or args.all_options:
        figdpi = args.sequence_logo_dpi
        figrowlength = args.sequence_logo_nt_per_row
        if args.sequence_logo_type == 'bar':
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)

            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=figdpi,
                                           figrowlength=figrowlength)
        elif args.sequence_logo_type == 'text':
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=figdpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=figrowlength)
        elif args.sequence_logo_type == 'both':
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)
            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=figdpi,
                                           figrowlength=figrowlength)
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=figdpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=figrowlength)

    if args.unalign_input:
        log.info("Generating a gap free version of the input alignment")
        if not args.silent:
            print("Generating a gap free version of the input alignment")
        outf = "%s_unaligned_input.fasta" % (args.outfile_stem)
        unaligned_arr = utilityFunctions.unAlign(orig_arr)
        utilityFunctions.writeOutfile(outf, unaligned_arr,
                                      orig_nams,
                                      removed_seqs)
    if args.replace_input:
        log.info("Generating a T instead of U version of the input alignment")
        if not args.silent:
            print("Generating a T instead of U version of the input alignment")
        outf = "%s_T_input.fasta" % (args.outfile_stem)
        T_arr = utilityFunctions.replaceUbyT(orig_arr)
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)

    if args.replace_output:
        log.info("Generating a T instead of U version of the output alignment")
        if not args.silent:
            print("Generating a T instead of U version of the output alignment")
        outf = "%s_T_output.fasta" % (args.outfile_stem)
        T_arr = utilityFunctions.replaceUbyT(arr)
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)

    if args.unalign_output:
        log.info("Generating a gap free version of the output alignment")
        if not args.silent:
            print("Generating a gap free version of the output alignment")
        outf = "%s_unaligned_output.fasta" % (args.outfile_stem)
        unaligned_arr = utilityFunctions.unAlign(arr)
        utilityFunctions.writeOutfile(outf, unaligned_arr,
                                      nams,
                                      removed_seqs)

    utilityFunctions.writeOutfile(outfile, arr, orig_nams,
                                  removed_seqs, rmfile)


if __name__ == "__main__":
    main()
