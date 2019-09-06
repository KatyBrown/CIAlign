#! /usr/bin/env python

import logging
import configargparse
import numpy as np
import copy
try:
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.parsingFunctions as parsingFunctions
    import CIAlign.miniAlignments as miniAlignments
    import CIAlign.similarityMatrix as similarityMatrix
    import CIAlign.consensusSeq as consensusSeq
except ModuleNotFoundError:
    import utilityFunctions
    import parsingFunctions
    import miniAlignments
    import similarityMatrix
    import consensusSeq


def main():
    parser = configargparse.ArgumentParser(
                description='''Clean and interpret a multiple sequence alignment''',
                add_help=False)
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')

    # Files
    # not to confuse with inifile
    required.add("--infile", dest='infile', type=str, required=True,
                 help='Path to input alignment. Required')
    optional.add("--inifile", dest='inifile', type=str,
                 default=None,
                 help='Path to config file. consensusSeq.pyDefault: %(default)s',
                 is_config_file=True)
    optional.add("--outfile_stem", dest='outfile_stem', type=str,
                 default="CIAlign",
                 help="Stem for output files (including path). Default: %(default)s")

    # Runtime
    optional.add("--silent", dest='silent',
                 help="Do not show progress on screen. Default: %(default)s",
                 action='store_true')

    # Crop Ends
    optional.add("--crop_ends", dest="crop_ends",
                 action="store_true",
                 help="Crop badly aligned ends. Default: %(default)s")
    optional.add("--crop_ends_mingap", dest='crop_ends_mingap',
                 type=int, default=10,
                 help="Minimum gap size to crop from ends. Default: %(default)s")

    # Remove Badly Aligned
    optional.add("--remove_badlyaligned", dest="remove_badlyaligned",
                 action="store_true",
                 help="Remove badly aligned sequences. Default: %(default)s")
    optional.add("--remove_badlyaligned_minperc", dest="remove_badlyaligned_minperc",
                 type=float, default=0.9,
                 help="Minimum percentage identity to majority to not be removed. Default: %(default)s")

    # Remove Insertions
    optional.add("--remove_insertions", dest="remove_insertions",
                 help="Remove insertions. Default: %(default)s",
                 action="store_true")
    optional.add("--insertion_min_size", dest="insertion_min_size",
                 type=int, default=3,
                 help="Minimum size insertion to remove. Default: %(default)s")
    optional.add("--insertion_max_size", dest="insertion_max_size",
                 type=int, default=300,
                 help="Maximum size insertion to remove. Default: %(default)s")
    optional.add("--insertion_min_flank", dest="insertion_min_flank",
                 type=int, default=5,
                 help="Minimum number of bases on either side of deleted insertions. Default: %(default)s")

    # Remove Short
    optional.add("--remove_short", dest="remove_short",
                 help="Remove sequences with less than n non-gap positions. Default: %(default)s",
                 action="store_true")
    optional.add("--remove_min_length", dest="remove_min_length",
                 type=int, default=50,
                 help="Minimum length sequence to remove. Default: %(default)s")

    optional.add("--remove_gaponly", dest="remove_gaponly",
                 action="store_false",
                 help="Remove gap only columns from the alignment. Default: %(default)s")

    # Consensus
    optional.add("--make_consensus", dest="make_consensus",
                 action="store_true",
                 help="Make a consensus sequence. Default: %(default)s")
    optional.add("--consensus_type", dest="consensus_type", type=str,
                 default="majority",
                 help="Type of consensus sequence to make. Default: %(default)s")
    optional.add("--consensus_keep_gaps", dest="consensus_keep_gaps",
                 action="store_true",
                 help="Keep gaps in consensus at positions where a gap is the consensus. Default: %(default)s")
    optional.add("--consensus_name", dest="consensus_name",
                 type=str, default="consensus",
                 help="Name of consensus sequence. Default: %(default)s")

    # Mini Alignments
    optional.add("--plot_input", dest="plot_input",
                 action="store_true",
                 help="Plot a mini alignment - an image representing the input alignment. Default: %(default)s")
    optional.add("--plot_output", dest="plot_output",
                 action="store_true",
                 help="Plot a mini alignment, an image reprsenting the output alignment. Default: %(default)s")
    optional.add("--plot_markup", dest="plot_markup",
                 action="store_true",
                 help="Plot the changes made to the alignment. Default: %(default)s")
    optional.add("--plot_dpi", dest="plot_dpi",
                 type=int, default=300,
                 help="DPI for mini alignments. Default: %(default)s")
    optional.add("--plot_format", dest="plot_format",
                 type=str, default='png',
                 help="Mini alignment plot format (png or svg). Default: %(default)s")
    optional.add("--plot_width", dest="plot_width",
                 type=int, default=5,
                 help="Width for mini alignments (inches). Default: %(default)s")
    optional.add("--plot_height", dest="plot_height",
                 type=int, default=3,
                 help="Height for mini alignments (inches). Default: %(default)s")

    # Sequence logos
    optional.add("--make_sequence_logo", dest="make_sequence_logo",
                 action="store_true",
                 help="Draw a sequence logo. Default: %(default)s")
    optional.add("--sequence_logo_type", dest="sequence_logo_type",
                 type=str, default='bar',
                 help="Type of sequence logo - bar/text/both. Default: %(default)s")
    optional.add("--sequence_logo_dpi", dest="sequence_logo_dpi",
                 type=int, default=300,
                 help="dpi for sequence logo image. Default: %(default)s")
    optional.add("--sequence_logo_font", dest="sequence_logo_font",
                 type=str, default='monospace',
                 help="Font for text sequence logo. Default: %(default)s")
    optional.add("--sequence_logo_nt_per_row", dest='sequence_logo_nt_per_row',
                 type=int, default=50,
                 help="Number of nucleotides or aas to show per row in the sequence logo. Default: %(default)s")
    optional.add("--sequence_logo_filetype", dest='sequence_logo_filetype',
                 type=str, default='png',
                 help="Output file type for sequence logo - png/jpg/svg. Default: %(default)s")
    optional.add("--list_fonts_only", dest='list_fonts_only',
                 action="store_true",
                 help="Make a swatch showing available fonts. Default: %(default)s")

    # Coverage
    optional.add("--plot_coverage_input", dest="plot_coverage_input",
                 action="store_true",
                 help="Plot the coverage of the input file as an interpolated function. Default: %(default)s")
    optional.add("--plot_coverage_output", dest="plot_coverage_output",
                 action="store_true",
                 help="Plot the coverage of the output file as an interpolated function. Default: %(default)s")
    optional.add("--plot_coverage_dpi", dest="plot_coverage_dpi",
                 type=int, default=300,
                 help="DPI for coverage plot. Default: %(default)s")
    optional.add("--plot_coverage_height", dest="plot_coverage_height",
                 type=int, default=3,
                 help="Height for coverage plot. Default: %(default)s")
    optional.add("--plot_coverage_width", dest="plot_coverage_width",
                 type=int, default=5,
                 help="Width for coverage plot. Default: %(default)s")
    optional.add("--plot_coverage_colour", dest="plot_coverage_colour",
                 type=str, default='#007bf5',
                 help="Colour for coverage plot. Default: %(default)s")
    optional.add("--plot_coverage_filetype", dest="plot_coverage_filetype",
                 type=str, default='png',
                 help="File type for coverage plot - can be png, jpg, tiff, svg. Default: %(default)s")

    # Similarity Matrix
    optional.add("--make_similarity_matrix_input", dest="make_simmatrix_input",
                 action="store_true",
                 help="Make a similarity matrix for the input alignment. Default: %(default)s")
    optional.add("--make_similarity_matrix_output", dest="make_simmatrix_output",
                 action="store_true",
                 help="Make a similarity matrix for the output alignment. Default: %(default)s")
    optional.add("--make_simmatrix_dp", dest="make_simmatrix_dp",
                 type=int, default=4,
                 help="N decimal places for the similarity matrix (output file only). Default: %(default)s")
    optional.add("--make_simmatrix_minoverlap",
                 dest="make_simmatrix_minoverlap",
                 type=int, default=1,
                 help="Minimum overlap between two sequences to have non-zero similarity in the similarity matrix. Default: %(default)s")
    optional.add("--make_simmatrix_keepgaps", dest="make_simmatrix_keepgaps",
                 type=bool, default=False,
                 help="Include positions with gaps in either or both sequences in the similarity matrix calculation. Default: %(default)s")

    optional.add('-h', '--help', action='help', default=configargparse.SUPPRESS,
                 help='Show this help message and exit')
    args = parser.parse_args()

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s_log.txt" % args.outfile_stem
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    rmfile = "%s_removed.txt" % args.outfile_stem
    outfile = "%s_parsed.fasta" % (args.outfile_stem)

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
        raise RuntimeError("Input alignment must be provided")

    arr, nams = utilityFunctions.FastaToArray(args.infile)

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
    outfile = "%s_parsed.fasta" % (args.outfile_stem)
    reset_rmfile = open(rmfile, "w")
    reset_rmfile.close()
    utilityFunctions.checkArrLength(arr, log)

    if args.crop_ends:
        # doesn't remove any whole columns or rows
        log.info("Cropping ends")
        if not args.silent:
            print("Cropping ends")
        arr, r = parsingFunctions.cropEnds(arr, nams, rmfile,
                                           log, args.crop_ends_mingap)
        markupdict['crop_ends'] = r
        removed_positions.update(r)
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_badlyaligned:
        log.info("Removing badly aligned sequences")
        if not args.silent:
            print("Removing badly aligned sequences")
        arr, r = parsingFunctions.removeBadlyAligned(arr, nams,
                                                     rmfile, log,
                                                     args.remove_badlyaligned_minperc)

        markupdict['remove_badlyaligned'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_insertions:
        log.info("Removing insertions")
        if not args.silent:
            print("Removing insertions")
        arr, r, relativePositions = parsingFunctions.removeInsertions(arr,
                                                                      relativePositions,
                                                                      rmfile,
                                                                      log,
                                                                      args.insertion_min_size,
                                                                      args.insertion_max_size,
                                                                      args.insertion_min_flank)

        markupdict['remove_insertions'] = r
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_short:
        log.info("Removing short sequences")
        if not args.silent:
            print("Removing short sequences")
        arr, r = parsingFunctions.removeTooShort(arr, nams, rmfile, log,
                                                 args.remove_min_length)

        markupdict['remove_short'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)
        utilityFunctions.checkArrLength(arr, log)

    if args.remove_gaponly:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        arr, r, relativePositions = parsingFunctions.removeGapOnly(arr,
                                                                   relativePositions,
                                                                   rmfile,
                                                                   log)
        markupdict['remove_gaponly'] = r
        utilityFunctions.checkArrLength(arr, log)

    if args.make_simmatrix_input:
        log.info("Building similarity matrix for input alignment")
        if not args.silent:
            print("Building similarity matrix for input alignment")
        outf = "%s_input_similarity.tsv" % (args.outfile_stem)
        similarityMatrix.calculateSimilarityMatrix(orig_arr,
                                                   orig_nams,
                                                   minoverlap=args.make_simmatrix_minoverlap,
                                                   keepgaps=args.make_simmatrix_keepgaps,
                                                   outfile=outf, dp=args.make_simmatrix_dp)

    if args.make_simmatrix_output:
        log.info("Building similarity matrix for output alignment")
        if not args.silent:
            print("Building similarity matrix for output alignment")
        outf = "%s_output_similarity.tsv" % (args.outfile_stem)
        similarityMatrix.calculateSimilarityMatrix(arr,
                                                   nams,
                                                   minoverlap=args.make_simmatrix_minoverlap,
                                                   keepgaps=args.make_simmatrix_keepgaps,
                                                   outfile=outf, dp=args.make_simmatrix_dp)

    if args.plot_input:
        log.info("Plotting mini alignment for input")
        if not args.silent:
            print("Building similarity matrix for output alignment")
        outf = "%s_input.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, orig_nams, log,
                                         outf, typ, args.plot_dpi,
                                         False, args.plot_width,
                                         args.plot_height)

    if args.plot_output:
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

    if args.plot_markup:
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

    if args.make_consensus:
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

    if args.plot_coverage_input:
        log.info("Plotting coverage for input")
        if not args.silent:
            print("Plotting coverage for input")
        coverage_file = "%s_input_coverage.%s" % (args.outfile_stem,
                                                  args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(orig_arr,
                                                     args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)

    if args.plot_coverage_output:
        log.info("Plotting coverage for output")
        if not args.silent:
            print("Plotting coverage for output")
        coverage_file = "%s_output_coverage.%s" % (args.outfile_stem,
                                                   args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(arr, args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)

    if args.make_sequence_logo:
        if args.sequence_logo_type == 'bar':
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
                out = "%s_logo_bar.%s" % (args.outfile_stem,
                                          args.sequence_logo_filetype)
                consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                               figdpi=args.sequence_logo_dpi,
                                               figrowlength=args.sequence_logo_nt_per_row)
        elif args.sequence_logo_type == 'text':
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=args.sequence_logo_dpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=args.sequence_logo_nt_per_row)
        elif args.sequence_logo_type == 'both':
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)
            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=args.sequence_logo_dpi,
                                           figrowlength=args.sequence_logo_nt_per_row)
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=args.sequence_logo_dpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=args.sequence_logo_nt_per_row)

    utilityFunctions.writeOutfile(outfile, arr, orig_nams,
                                  removed_seqs, rmfile)

if __name__ == "__main__":
    main()
