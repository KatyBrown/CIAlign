#! /usr/bin/env python

import logging
import configargparse
import os.path
import numpy as np
import copy
import argP

try:
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.parsingFunctions as parsingFunctions
    import CIAlign.miniAlignments as miniAlignments
    import CIAlign.similarityMatrix as similarityMatrix
    import CIAlign.consensusSeq as consensusSeq
except ImportError:
    import utilityFunctions
    import parsingFunctions
    import miniAlignments
    import similarityMatrix
    import consensusSeq


def main():
    # Get and parse the argument parser
    parser = argP.getParser()
    args = parser.parse_args()

    # Set up logger
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s_log.txt" % args.outfile_stem
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # Create a logging format
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

    arr, nams = utilityFunctions.FastaToArray(args.infile, log,
                                              args.outfile_stem)
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
        assert args.insertion_min_size < args.insertion_max_size, "\
            insertion_min_size must be less than insertion_max_size"
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
    fn = args.plot_force_numbers
    if args.plot_input or args.all_options:
        log.info("Plotting mini alignment for input")
        if not args.silent:
            print("Plotting mini alignment for input")
        outf = "%s_input.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, log, cmt_file,
                                         outf, typ, args.plot_dpi,
                                         False, args.plot_width,
                                         args.plot_height,
                                         force_numbers=fn)

    if args.plot_output or args.all_options:
        log.info("Plotting mini alignment for output")
        if not args.silent:
            print("Plotting mini alignment for output")
        outf = "%s_output.%s" % (args.outfile_stem, args.plot_format)
        if not args.plot_keep_numbers:
            miniAlignments.drawMiniAlignment(arr, nams, log,
                                             outf, typ,
                                             args.plot_dpi,
                                             False,
                                             args.plot_width,
                                             args.plot_height,
                                             force_numbers=fn)
        else:
            miniAlignments.drawMiniAlignment(arr, nams, log,
                                             outf, typ,
                                             args.plot_dpi,
                                             False,
                                             args.plot_width,
                                             args.plot_height,
                                             orig_nams=orig_nams,
                                             keep_numbers=True,
                                             force_numbers=fn)

    if args.plot_markup or args.all_options:
        log.info("Plotting mini alignment with markup")
        if not args.silent:
            print("Plotting mini alignment with markup")
        outf = "%s_markup.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, orig_nams, log,
                                         outf, typ,
                                         args.plot_dpi,
                                         False,
                                         args.plot_width,
                                         args.plot_height,
                                         markup=True,
                                         markupdict=markupdict,
                                         force_numbers=fn)

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
        log.info("Generating a T instead of U version of\
                 the output alignment")
        if not args.silent:
            print("Generating a T instead of U version of\
                  the output alignment")
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
