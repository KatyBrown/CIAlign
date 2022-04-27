#!/usr/bin/env python3
import os
import copy
import numpy as np
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


def run(args, log):

    # Basic checks before running
    prelimChecks(args, log)
    # Set up arrays of the sequence names and aligned sequences
    arr, nams, typ = setupArrays(args, log)

    # Make copies of the unedited arrays
    orig_arr = copy.copy(arr)
    orig_nams = copy.copy(nams)
    functions = whichFunctions(args)

    if "cleaning" in functions:
        arr, nams, markupdict, removed = runCleaning(args,
                                                     log,
                                                     arr,
                                                     nams)
    else:
        markupdict = dict()
        removed = set()

    if "matrices" in functions:
        # Make similarity matrices
        runMatrix(args, log, orig_arr, orig_nams, arr, nams)

    if "mini_alignments" in functions:
        # Plot mini alignments
        runMiniAlignments(args, log, orig_arr, orig_nams, arr, nams,
                          markupdict, typ)

    if "consensus" in functions:
        # Make consensus sequences
        runConsensus(args, log, orig_arr, orig_nams, arr, nams, removed)

    if "coverage" in functions:
        # Plot coverage plots
        runCoverage(args, log, orig_arr, orig_nams, arr, nams)

    if "logos" in functions:
        # Draw sequence logos
        runSeqLogo(args, log, orig_arr, orig_nams, arr, nams, typ)

    if "unalign" in functions:
        # Make an unaligned copy of the input or output
        runUnalign(args, log, orig_arr, orig_nams, arr, nams, removed)

    if "ttou" in functions:
        # Convert T to U in the input or output
        runTtoU(args, log, orig_arr, orig_nams, arr, nams, removed)


def prelimChecks(args, log):
    '''
    Initial checks that everything is OK - run list_fonts_only of that is
    specified, otherwise check the input alignment can be found.
    Exits if anything is amiss.

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    '''
    # If list_fonts_only is specified CIAlign will just make an image with
    # all available fonts and exit
    if args.list_fonts_only:
        print("Listing fonts")
        out = "%s_fonts.png" % args.outfile_stem
        utilityFunctions.listFonts(out)
        exit()

    # Check the input file exists and can be found

    # XXXXX comes from the ini file
    if not args.infile or args.infile == 'XXXXX':
        print("Input alignment must be provided.")
        exit()

    # check if input file exists and is a file
    if not os.path.isfile(args.infile):
        print("Error! Your input alignmnent path could not be found.")
        exit()


def whichFunctions(args):
    '''
    Make a list to track which groups of functions to run

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters

    Returns
    -------
    which_functions: list
        List of groups of functions to run
    '''
    which_functions = []
    # Cleaning Functions
    if any([args.remove_divergent,
            args.remove_insertions,
            args.crop_ends,
            args.remove_short,
            args.remove_gaponly,
            args.clean,
            args.all_options]):
        which_functions.append("cleaning")

    # Similarity Matrix
    if any([args.make_simmatrix_input,
            args.make_simmatrix_output,
            args.interpret,
            args.all_options]):
        which_functions.append("matrices")

    # Mini Alignments
    if any([args.plot_input,
            args.plot_output,
            args.plot_markup,
            args.visualise,
            args.all_options]):
        which_functions.append("mini_alignments")

    # Consensus sequences
    if any([args.make_consensus,
            args.interpret,
            args.all_options]):
        which_functions.append("consensus")

    # Coverage Plots
    if any([args.plot_coverage_input,
            args.plot_coverage_output,
            args.interpret,
            args.all_options]):
        which_functions.append("coverage")

    # Sequence Logos
    if any([args.make_sequence_logo,
            args.all_options]):
        which_functions.append("logos")

    # Unalign input or output
    if any([args.unalign_input,
            args.unalign_output,
            args.all_options]):
        which_functions.append("unalign")

    # Replace T with U in input or output
    if any([args.replace_input,
            args.replace_output]):
        which_functions.append("ttou")

    return (which_functions)


def setupArrays(args, log):
    '''
    Read the alignment into an array, check there are enough sequences in
    the array and the names are not duplicated, detect if the alignment
    is nucleotides or amino acids.

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file

    Returns
    -------
    arr: np.array
        The alignment stored in a numpy array
    nams: list
        The names of the sequences in the alignment
    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide
    '''

    # convert the input fasta file into an array and make a list of
    # sequence names so the order can be maintained
    arr, nams = utilityFunctions.FastaToArray(args.infile, log,
                                              args.outfile_stem)
    # check if names are unique
    if len(nams) > len(set(nams)):
        print("Error! Your input alignmnent has duplicate names!")
        exit()

    # Check the alignment array isn't empty
    utilityFunctions.checkArrLength(arr, log)

    # Check which cleaning functions are requested
    cleaningArgs = [args.remove_insertions,
                    args.crop_ends,
                    args.remove_divergent]

    # Check there are enough sequences for the requested functions
    if len(arr) < 3 and any(cleaningArgs):
        # when less than three sequences, stop
        print("You need at least three sequences in your MSA to run \
               remove_insertions, crop_ends or remove_divergent")
        exit()
    elif len(arr) < 2:
        print("You need at least two sequences in your MSA")
        exit()

    # detect if the sequence is amino acids or nucleotides
    typ = utilityFunctions.seqType(arr)

    if typ == 'aa':
        log.info("Amino acid alignment detected")
    else:
        log.info("Nucleotide alignment detected")

    return (arr, nams, typ)


def setupTrackers(args, arr):
    '''
    Sets up variables to store the rows, columns and postions removed
    by the cleaning functions.

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    arr: np.array
        The alignment stored in a numpy array

    Returns
    -------
    markupdict: dict
        Dictionary where the keys will be function names and the values are
        lists of columns, rows or positions which have been removed
    relativePositions: list
        A list of integers representing columns in the alignment, from which
        values are removed as columns are removed from the alignment.
    removed_seqs: set
        set to store the names of sequences which have been removed
    removed_cols: set
        set to store the indices of columns which have been removed
    removed_positions: dict
        Dictionary which will have sequence names as keys and tuples as values,
        where tuple[0] is a list of positions which have been removed at the
        beginning of the sequence and tuple[1] is a list of positions which
        have been removed at the end of the sequence
    '''
    # make a dictionary to store the changes made
    markupdict = dict()

    # remember positions relative to original alignment
    relativePositions = list(range(0, len(arr[0])))

    removed_seqs = set()
    removed_cols = set()
    removed_positions = dict()

    return (markupdict, relativePositions,
            [removed_seqs, removed_cols, removed_positions])


def setupOutfiles(args):
    '''
    Set up output files for the cleaning functions

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters

    Returns
    -------
    outfile: str
        Path to the main output file for the cleaned fasta file
    rmfile: str
        Path to a file to store the removed columns, rows and positions
        formatted to be parsed automatically.
    '''
    # File for removed sequence list
    rmfile = "%s_removed.txt" % args.outfile_stem
    # File for cleaned alignment
    outfile = "%s_cleaned.fasta" % (args.outfile_stem)
    # Remove any text which is already in the rmfile
    reset_rmfile = open(rmfile, "w")
    reset_rmfile.close()
    return (outfile, rmfile)


def runCleaning(args, log, arr, nams):
    '''
    Run the cleaning functions

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    arr: np.array
        Array containing the original alignment
    nams:
        List of sequence names in the original alignment

    Returns
    -------
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names remaining in the cleaned alignment
    markupdict: dict
        Dictionary where the keys are function names and the values are
        lists of columns, rows or positions which have been removed
    removed_seqs: set
        set of the names of sequences which have been removed
    '''
    # Set everything up
    orig_nams = copy.copy(nams)
    markupdict, relativePositions, R = setupTrackers(args, arr)
    outfile, rmfile = setupOutfiles(args)
    removed_seqs, removed_cols, removed_positions = R

    # Remove divergent sequences
    if args.remove_divergent or args.all_options or args.clean:
        log.info("Removing divergent sequences")
        if not args.silent:
            print("Removing divergent sequences")
        minperc = args.remove_divergent_minperc
        arr, r = parsingFunctions.removeDivergent(arr, nams,
                                                  rmfile, log,
                                                  minperc)
        # Track what has been removed
        markupdict['remove_divergent'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)

        # Check there are some sequences left
        utilityFunctions.checkArrLength(arr, log)

    # Remove gaps created by remove divergent
    if (args.remove_divergent
            and args.remove_gaponly) or args.all_options or args.clean:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")
        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        # Track what has been removed
        arr, r, relativePositions = A

        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r

        # Check there are some columns left
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    # Remove insertions
    if args.remove_insertions or args.all_options or args.clean:
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
                                              args.insertion_min_flank,
                                              args.insertion_min_perc)

        # Track what has been removed
        arr, r, relativePositions = A
        markupdict['remove_insertions'] = r
        removed_cols = removed_cols | r
        # Check there are some columns left
        utilityFunctions.checkArrLength(arr, log)

    # Remove gaps created by remove insertions
    if (args.remove_insertions
            and args.remove_gaponly) or args.all_options or args.clean:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")
        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)

        # Track what has been removed
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        # Check there are still some columns left
        utilityFunctions.checkArrLength(arr, log)

    # Crop Ends
    if args.crop_ends or args.all_options or args.clean:
        # doesn't remove any whole columns or rows
        log.info("Cropping ends")
        if not args.silent:
            print("Cropping ends")
        arr, r = parsingFunctions.cropEnds(arr, nams, relativePositions,
                                           rmfile,
                                           log, args.crop_ends_mingap_perc,
                                           args.crop_ends_redefine_perc)
        # Track what has been removed
        markupdict['crop_ends'] = r
        removed_positions.update(r)
        # Check there are still some positions left
        utilityFunctions.checkArrLength(arr, log)

    # Remove empty columns created by crop ends
    if (args.crop_ends
            and args.remove_gaponly) or args.all_options or args.clean:
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        # Track what has been removed
        arr, r, relativePositions = A
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        # Check there are still some positions left
        utilityFunctions.checkArrLength(arr, log)

    # Remove short
    if args.remove_short or args.all_options or args.clean:
        log.info("Removing short sequences")
        if not args.silent:
            print("Removing short sequences")
        arr, r = parsingFunctions.removeTooShort(arr, nams, rmfile, log,
                                                 args.remove_min_length)
        # Track what has been removed
        markupdict['remove_short'] = r
        removed_seqs = removed_seqs | r
        nams = utilityFunctions.updateNams(nams, r)
        # Check there are still some sequences left
        utilityFunctions.checkArrLength(arr, log)

    # Remove empty columns created by remove short
    if (args.remove_short
            and args.remove_gaponly) or args.all_options or args.clean:
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
                                    args.remove_short or
                                    args.clean):
        log.info("Removing gap only columns")
        if not args.silent:
            print("Removing gap only columns")

        A = parsingFunctions.removeGapOnly(arr,
                                           relativePositions,
                                           rmfile,
                                           log)
        arr, r, relativePositions = A
        # Track what has been removed
        if 'remove_gaponly' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r
        removed_cols = removed_cols | r
        # Check there are some columns left
        utilityFunctions.checkArrLength(arr, log)

    # Write the output file
    utilityFunctions.writeOutfile(outfile, arr, orig_nams,
                                  removed_seqs, rmfile)

    return (arr, nams, markupdict, removed_seqs)


def runMatrix(args, log, orig_arr, orig_nams, arr, nams):
    '''
    Make similarity matrices

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams:
        List of sequence names in the cleaned alignment
    '''
    # Input matrix
    if args.make_simmatrix_input or args.all_options or args.interpret:
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
    # Output matrix
    # todo: what if only interpret functions are called?
    if args.make_simmatrix_output or args.all_options or args.interpret:
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


def runMiniAlignments(args, log, orig_arr, orig_nams, arr, nams,
                      markupdict, typ):
    '''
    Plot mini alignments

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    markupdict: dict
        Dictionary where the keys are function names and the values are
        lists of columns, rows or positions which have been removed
    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide
    '''
    fn = args.plot_force_numbers
    # Mini alignment of CIAlign input
    if args.plot_input or args.all_options or args.visualise:
        log.info("Plotting mini alignment for input")
        if not args.silent:
            print("Plotting mini alignment for input")
        outf = "%s_input.%s" % (args.outfile_stem, args.plot_format)
        miniAlignments.drawMiniAlignment(orig_arr, orig_nams, log,
                                         outf, typ, args.plot_dpi,
                                         False, args.plot_width,
                                         args.plot_height,
                                         force_numbers=fn)
    # Mini alignment of CIAlign output
    # todo: what if only interpret functions are called?
    if args.plot_output or args.all_options or args.visualise:
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
    # Markup plot
    # todo: what if only interpret functions are called?
    if args.plot_markup or args.all_options or args.visualise:
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


def runConsensus(args, log, orig_arr, orig_nams, arr, nams, removed_seqs):
    '''
    Make consensus sequences

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    removed_seqs: set
        Set of sequence names which have been removed
    '''
    if args.make_consensus or args.all_options or args.interpret:
        log.info("Building consensus sequence")
        if not args.silent:
            print("Building consensus sequence")
        cons, coverage = consensusSeq.findConsensus(arr,
                                                    log, args.consensus_type)
        consarr = np.array(cons)
        # Combine the consensus with the alignment
        arr_plus_cons = np.row_stack((arr, consarr))
        cons = "".join(cons)
        # Remove the gaps from the consensus if this option is specified
        if not args.consensus_keep_gaps:
            cons = cons.replace("-", "")
        # Output file of just the consensus sequence
        out = open("%s_consensus.fasta" % args.outfile_stem, "w")
        out.write(">%s\n%s\n" % (args.consensus_name, cons))
        out.close()
        # Output file of the consensus and the alignment
        outf = "%s_with_consensus.fasta" % args.outfile_stem
        utilityFunctions.writeOutfile(outf, arr_plus_cons,
                                      nams + [args.consensus_name],
                                      removed_seqs)


def runCoverage(args, log, orig_arr, orig_nams, arr, nams):
    '''
    Draw coverage plots

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    '''
    # Coverage plot for CIAlign input
    if args.plot_coverage_input or args.all_options or args.interpret:
        log.info("Plotting coverage for input")
        if not args.silent:
            print("Plotting coverage for input")
        coverage_file = "%s_input_coverage.%s" % (args.outfile_stem,
                                                  args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(orig_arr,
                                                     args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)

    # Coverage plot for CIAlign output
    # todo: what if only interpret functions are called?
    if args.plot_coverage_output or args.all_options or args.interpret:
        if not args.silent:
            print("Plotting coverage for output")
        log.info("Plotting coverage for output")
        coverage_file = "%s_output_coverage.%s" % (args.outfile_stem,
                                                   args.plot_coverage_filetype)
        consx, coverage = consensusSeq.findConsensus(arr, args.consensus_type)
        consensusSeq.makeCoveragePlot(coverage, coverage_file)


def runSeqLogo(args, log, orig_arr, orig_nams, arr, nams, typ):
    '''
    Plot sequence logos

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide
    '''
    if args.make_sequence_logo or args.all_options:
        figdpi = args.sequence_logo_dpi
        figrowlength = args.sequence_logo_nt_per_row
        logo_start = args.logo_start
        logo_end = args.logo_end
        if logo_end < logo_start:
            print("Error! The start should be smaller than the end for the \
                  sequence logo!")
            exit()
        # Sequence logo bar chart
        if args.sequence_logo_type == 'bar':
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)
            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=figdpi,
                                           figrowlength=figrowlength,
                                           start=logo_start, end=logo_end)
        elif args.sequence_logo_type == 'text':
            # Text sequence logo
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=figdpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=figrowlength,
                                       start=logo_start, end=logo_end)
        elif args.sequence_logo_type == 'both':
            # Plot both types of sequence logo
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)
            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=figdpi,
                                           figrowlength=figrowlength,
                                           start=logo_start, end=logo_end)
            log.info("Generating text sequence logo")
            if not args.silent:
                print("Generating text sequence logo")
            out = "%s_logo_text.%s" % (args.outfile_stem,
                                       args.sequence_logo_filetype)
            consensusSeq.sequence_logo(arr, out, typ=typ,
                                       figdpi=figdpi,
                                       figfontname=args.sequence_logo_font,
                                       figrowlength=figrowlength,
                                       start=logo_start, end=logo_end)


def runUnalign(args, log, orig_arr, orig_nams, arr, nams, removed_seqs):
    '''
    Make a copy of the alignment without gaps

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    removed_seqs: set
        Set of sequence names which have been removed
    '''
    # Unalign input
    if args.unalign_input:
        log.info("Generating a gap free version of the input alignment")
        if not args.silent:
            print("Generating a gap free version of the input alignment")
        outf = "%s_unaligned_input.fasta" % (args.outfile_stem)
        unaligned_arr = utilityFunctions.unAlign(orig_arr)
        # Write to file
        utilityFunctions.writeOutfile(outf, unaligned_arr,
                                      orig_nams,
                                      removed_seqs)
    # Unalign output
    if args.unalign_output:
        log.info("Generating a gap free version of the output alignment")
        if not args.silent:
            print("Generating a gap free version of the output alignment")
        outf = "%s_unaligned_output.fasta" % (args.outfile_stem)
        unaligned_arr = utilityFunctions.unAlign(arr)
        # Write to file
        utilityFunctions.writeOutfile(outf, unaligned_arr,
                                      nams,
                                      removed_seqs)


def runTtoU(args, log, orig_arr, orig_nams, arr, nams, removed_seqs):
    '''
    Make a copy of the alignment with T replaced by U

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    orig_arr: np.array
        Array containing the original alignment
    orig_nams:
        List of sequence names in the original alignment
    arr: np.array
        Array containing the cleaned alignment
    nams: list
        List of sequence names in the cleaned alignment
    removed_seqs: set
        Set of sequence names which have been removed
    '''
    # Replace T with U in the input
    if args.replace_input:
        log.info("Generating a T instead of U version of the input alignment")
        if not args.silent:
            print("Generating a T instead of U version of the input alignment")
        outf = "%s_T_input.fasta" % (args.outfile_stem)
        T_arr = utilityFunctions.replaceUbyT(orig_arr)
        # Write to file
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)
    # Rpleace T with U in the output
    if args.replace_output:
        log.info("Generating a T instead of U version of\
                 the output alignment")
        if not args.silent:
            print("Generating a T instead of U version of\
                  the output alignment")
        outf = "%s_T_output.fasta" % (args.outfile_stem)
        T_arr = utilityFunctions.replaceUbyT(arr)
        # Write to file
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)
