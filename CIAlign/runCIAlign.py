#!/usr/bin/env python3
import os
import copy
import numpy as np
import pandas as pd
try:
    import CIAlign.utilityFunctions as utilityFunctions
    import CIAlign.parsingFunctions as parsingFunctions
    import CIAlign.miniAlignments as miniAlignments
    import CIAlign.similarityMatrix as similarityMatrix
    import CIAlign.consensusSeq as consensusSeq
    import CIAlign.matrices as matrices
except ImportError:
    import utilityFunctions
    import parsingFunctions
    import miniAlignments
    import similarityMatrix
    import consensusSeq
    import matrices


def run(args, log):

    # Basic checks before running
    prelimChecks(args, log)
    # Set up arrays of the sequence names and aligned sequences
    arr, nams, typ = setupArrays(args, log)

    # Make copies of the unedited arrays
    orig_arr = copy.copy(arr)
    orig_nams = copy.copy(nams)
    functions = whichFunctions(args)

    if "section" in functions:
        arr, removed_c = setupSection(args, log, arr)
    else:
        removed_c = set()

    if "cleaning" in functions:
        keeps = setupRetains(args, nams, log)
        arr, nams, markupdict, removed_r, removed_c = runCleaning(args,
                                                                  log,
                                                                  orig_arr,
                                                                  arr,
                                                                  nams,
                                                                  keeps,
                                                                  removed_c)
    else:
        markupdict = dict()
        removed_c = set()
        removed_r = set()

    if "matrices" in functions:
        # Make similarity matrices
        runMatrix(args, log, orig_arr, orig_nams, arr, nams)

    if "mini_alignments" in functions:
        # Plot mini alignments
        runMiniAlignments(args, log, orig_arr, orig_nams, arr, nams,
                          markupdict, typ)

    if "consensus" in functions:
        # Make consensus sequences
        runConsensus(args, log, orig_arr, orig_nams, arr, nams, removed_r)

    if "stats" in functions:
        # Plot coverage plots
        runStatsPlots(args, log, orig_arr, orig_nams, arr, nams, typ)

    if "logos" in functions:
        # Draw sequence logos
        runSeqLogo(args, log, orig_arr, orig_nams, arr, nams, typ, removed_c)

    if "unalign" in functions:
        # Make an unaligned copy of the input or output
        runUnalign(args, log, orig_arr, orig_nams, arr, nams, removed_r)

    if "ttou" in functions:
        # Convert T to U in the input or output
        runTtoU(args, log, orig_arr, orig_nams, arr, nams, removed_r,
                rev=False)

    if "utot" in functions:
        # Convert U to T in the input or output
        runTtoU(args, log, orig_arr, orig_nams, arr, nams, removed_r,
                rev=True)

    if "pwm" in functions:
        # Generate position matrices
        runPWM(args, log, orig_arr, arr, typ, removed_c)


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
    if args.get_section:
        if args.section_start is None or args.section_end is None:
            print("Error! Start (--section_start) and end \
(--section_end) positions must be provided with --get_section")
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
            args.crop_divergent,
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
    if any([args.plot_stats_input,
            args.plot_stats_output,
            args.interpret,
            args.all_options]):
        which_functions.append("stats")

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
    if any([args.replace_input_tu,
            args.replace_output_tu]):
        which_functions.append("ttou")

    # Replace T with U in input or output
    if any([args.replace_input_ut,
            args.replace_output_ut,]):
        which_functions.append("utot")

    if any([args.pwm_input,
            args.pwm_output]):
        which_functions.append("pwm")

    if args.get_section:
        which_functions.append("section")
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
    # check if names are unique.append(
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
    typ = utilityFunctions.seqType(arr, log)

    if typ == 'aa':
        log.info("Amino acid alignment detected")
    else:
        log.info("Nucleotide alignment detected")

    return (arr, nams, typ)


def setupRetains(args, nams, log):
    '''
    Sets up a dictionary of sequence names which the cleaning functions should
    ignore.

    These come from the "retain" command line arguments and can currently
    be specified for crop_ends, remove_divergent, remove_short or all
    rowwise functions.

    Allows the user to specify sequences to keep regardless of whether
    they pass or fail the rowwise cleaning operation thresholds.

    Sequence names can be specified individually on the command line with
    --retain, --crop_ends_retain, --remove_divergent_retain,
    --remove_short_retain

    They can also be specified as a list in a text file with --retain_list,
    --crop_ends_retain_list etc, in which case the value is the path to
    the file.

    Finally they can be specified by searching each name for a character
    string, specified as --retain_str, --crop_ends_retain_str etc.

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    nams: list
        The names of the sequences in the alignment
    log: logging.Logger
        Open log file

    Returns
    -------
    keepD: dict
        A dictionary listing sequences not to process for each function,
        where the keys are function names.
    '''
    # These are all the arguments from the command line to retain
    # sequences, no suffix for directly specified, S for string, L for list,
    # rs = remove_short, rd = remove_divergent, ce = crop_ends
    retain_args = [args.retain_seqs_rs,
                   args.retain_seqs_rsS,
                   args.retain_seqs_rsL,
                   args.retain_seqs_rd,
                   args.retain_seqs_rdS,
                   args.retain_seqs_rdL,
                   args.retain_seqs_ce,
                   args.retain_seqs_ceS,
                   args.retain_seqs_ceL,
                   args.retain_seqs,
                   args.retain_seqsS,
                   args.retain_seqsL]
    # Turn the list into an array
    retain_args = np.array(retain_args, dtype=object)
    keepD = dict()

    # Split the parameters into four batches of 3 variables
    rr = np.split(retain_args, 4)

    # These are the functions these parameters are used with, each corresponds
    # to a batch in rr
    titles = ['remove_short',
              'remove_divergent',
              'crop_ends',
              'all_rowwise']

    for i, r in enumerate(rr):
        # Make an array listing the sequences referenced by each variable
        keeps = utilityFunctions.configRetainSeqs(r[0],
                                                  r[1],
                                                  r[2],
                                                  nams,
                                                  titles[i],
                                                  log,
                                                  args.silent)
        # Store the resulting array in a dictionary where the key
        # is the function name
        keepD[titles[i]] = keeps
    return (keepD)


def setupTrackers(args, arr, removed_cols, rmfile):
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
    for col in removed_cols:
        relativePositions.remove(col)
    if len(removed_cols) != 0:
        markupdict['user'] = removed_cols
    removed_positions = dict()
    if len(removed_cols) != 0:
        out = open(rmfile, "a")
        out.write("user_defined\t%s\n" % (",".join(
            [str(r) for r in removed_cols])))
        out.close()
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


def setupSection(args, log, arr):
    '''
    Extracts a subset of columns from the alignment for processing, as
    specified by the user.

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        Open log file
    arr: np.array
        Array containing the original alignment

    Returns
    -------
    arr: np.array
        Array containing only the specified section of the alignment
    removed_c: set
        Set containing the indices of the removed columns
    '''
    orig_width = np.shape(arr)[1]
    assert args.section_end - args.section_start > 5, (
        "Section must be at least 5 residues in length")
    log.info("Cropping alignment to keep columns %i to %i" % (
        args.section_start, args.section_end))
    if not args.silent:
        print("Cropping alignment to keep columns %i to %i" % (
            args.section_start, args.section_end))

    arr = copy.copy(arr[:, args.section_start:args.section_end+1])
    removed_c = list(np.arange(0, args.section_start))
    removed_c += list(np.arange(args.section_end+1, orig_width))
    removed_c = set(removed_c)
    return (arr, removed_c)


def runCleaning(args, log, orig_arr, arr, nams, keeps, removed_c):
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
    outfile, rmfile = setupOutfiles(args)
    markupdict, relativePositions, R = setupTrackers(args, orig_arr,
                                                     removed_c, rmfile)

    removed_seqs, removed_cols, removed_positions = R

    # Remove divergent sequences
    if args.remove_divergent or args.all_options or args.clean:
        log.info("Removing divergent sequences")
        if not args.silent:
            print("Removing divergent sequences")
        minperc = args.remove_divergent_minperc
        arr, r = parsingFunctions.removeDivergent(arr, nams,
                                                  rmfile, log,
                                                  keeps,
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

        if 'remove_gap_only' in markupdict:
            markupdict['remove_gaponly'].update(r)
        else:
            markupdict['remove_gaponly'] = r

        # Check there are some columns left
        removed_cols = removed_cols | r
        utilityFunctions.checkArrLength(arr, log)

    # Crop divergent
    if args.crop_divergent:
        log.info("Removing divergent sequence ends")
        if not args.silent:
            print("Removing divergent sequence ends")

        A = parsingFunctions.cropDivergent(arr,
                                           relativePositions,
                                           rmfile,
                                           log,
                                           args.divergent_min_prop_ident,
                                           args.divergent_min_prop_nongap,
                                           args.divergent_buffer_size)

        # Track what has been removed
        arr, r, relativePositions = A
        markupdict['crop_divergent'] = r
        removed_cols = removed_cols | r
        # Check there are some columns left
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
        if 'remove_gap_only' in markupdict:
            markupdict['remove_gap_only'].update(r)
        else:
            markupdict['remove_gap_only'] = r
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
                                           log, keeps,
                                           args.crop_ends_mingap_perc,
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
                                                 keeps,
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
                                    args.crop_divergent or
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

    return (arr, nams, markupdict, removed_seqs, removed_cols)


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
                                         force_numbers=fn,
                                         palette=args.palette)
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
                                             force_numbers=fn,
                                             palette=args.palette)
        else:
            miniAlignments.drawMiniAlignment(arr, nams, log,
                                             outf, typ,
                                             args.plot_dpi,
                                             False,
                                             args.plot_width,
                                             args.plot_height,
                                             orig_nams=orig_nams,
                                             keep_numbers=True,
                                             force_numbers=fn,
                                             palette=args.palette)
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
                                         force_numbers=fn,
                                         palette=args.palette)


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


def runStatsPlots(args, log, orig_arr, orig_nams, arr, nams, typ):
    '''
    Draw plots of different statistics about the alignment.

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
    to_plot = []
    if args.plot_stats_input or args.all_options or args.interpret:
        to_plot.append("input")
    if args.plot_stats_output or args.all_options or args.interpret:
        to_plot.append("output")

    for inout in to_plot:
        if inout == 'input':
            c_arr = copy.deepcopy(orig_arr)
        else:
            c_arr = copy.deepcopy(arr)
        rowD = dict()
        consx, coverage = consensusSeq.findConsensus(c_arr,
                                                     args.consensus_type)
        rowD['coverage'] = coverage
        bit_scores, ents = consensusSeq.calcConservationAli(c_arr, typ)
        rowD['information_content'] = bit_scores
        rowD['shannon_entropy'] = ents

        log.info("Plotting coverage for %s" % inout)
        if not args.silent:
            print("Plotting coverage for %s" % inout)
        for stat in rowD:
            outfile = "%s_%s_%s.%s" % (args.outfile_stem, inout,
                                       stat, args.plot_stats_filetype)

            consensusSeq.makeLinePlot(rowD[stat],
                                      outfile,
                                      stat.replace("_", " ").title(),
                                      dpi=args.plot_stats_dpi,
                                      height=args.plot_stats_height,
                                      width=args.plot_stats_width,
                                      colour=args.plot_stats_colour)

        stats_tab = pd.DataFrame(rowD.values()).T.round(4)
        stats_tab.columns = rowD.keys()
        stats_tab.to_csv("%s_%s_column_stats.tsv" % (args.outfile_stem,
                                                     inout), sep="\t")


def runSeqLogo(args, log, orig_arr, orig_nams, arr, nams, typ, removed):
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
    if args.make_sequence_logo or args.all_options or args.visualise:
        figdpi = args.sequence_logo_dpi
        figrowlength = args.sequence_logo_nt_per_row
        logo_start, logo_end = utilityFunctions.updateStartEnd(
            args.logo_start, args.logo_end, removed)
        if logo_end < logo_start:
            print("Error! The start should be smaller than the end for the \
                  sequence logo!")
            exit()
        # Sequence logo bar chart
        if (args.sequence_logo_type == 'bar' or
                args.sequence_logo_type == 'both'):
            log.info("Generating sequence logo bar chart")
            if not args.silent:
                print("Generating sequence logo bar chart")
            out = "%s_logo_bar.%s" % (args.outfile_stem,
                                      args.sequence_logo_filetype)
            consensusSeq.sequence_bar_logo(arr, out, typ=typ,
                                           figdpi=figdpi,
                                           figrowlength=figrowlength,
                                           start=logo_start, end=logo_end,
                                           palette=args.palette)
        if (args.sequence_logo_type == 'text' or
                args.sequence_logo_type == 'both'):
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
                                       start=logo_start, end=logo_end,
                                       palette=args.palette)


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


def runTtoU(args, log, orig_arr, orig_nams, arr, nams, removed_seqs,
            rev=False):
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
    rev: bool
        If True, do the opposite, change U to T
    '''
    if not rev:
        one = "T"
        two = "U"
    else:
        one = "U"
        two = "T"

    # Replace T with U in the input
    if args.replace_input_tu or args.replace_input_ut:
        log.info("Generating a %s instead of %s version of the input \
alignment" % (two, one))
        if not args.silent:
            print("Generating a %s instead of %s version of the input \
alignment" % (two, one))
        outf = "%s_%s_input.fasta" % (args.outfile_stem, two)
        T_arr = utilityFunctions.replaceUbyT(orig_arr, rev=rev)
        # Write to file
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)
    # Replace T with U in the output
    if args.replace_output_tu or args.replace_output_ut:
        log.info("Generating a %s instead of %s version of the output \
alignment" % (two, one))
        if not args.silent:
            print("Generating a %s instead of %s version of the output \
alignment" % (two, one))
        outf = "%s_%s_output.fasta" % (args.outfile_stem, two)
        T_arr = utilityFunctions.replaceUbyT(arr, rev=rev)
        # Write to file
        utilityFunctions.writeOutfile(outf, T_arr,
                                      orig_nams,
                                      removed_seqs)


def runPWM(args, log, orig_arr, arr, typ, removed):
    '''
    Converts an alignment array into a PFM, PPM and PWM.

    PFM - position frequency matrix - a matrix showing the number of
    each residue in each column of the alignment.
    PPM - position probability matrix - normalised by background
    frequency and including pseudocounts.
    PWM - position weight matrix - a matrix which shows the log-likelihood
    ratio of observing character i at position j in a site
    compared with a random sequence (from 10.1186/s12859-020-3348-6)

    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    log: logging.Logger
        An open log file object
    arr: np.array
        The alignment the PFM is for, stored as a numpy array
    typ: str
        nt for nucleotide or aa for amino acid
    '''
    runs = []
    if args.pwm_input:
        runs.append('input')
    if args.pwm_output:
        runs.append('output')
    for run in runs:
        if run == 'input':
            c_arr = copy.deepcopy(orig_arr)
        elif run == 'output':
            c_arr = copy.deepcopy(arr)
        if args.pwm_start is not None and args.pwm_end is not None:
            start = args.pwm_start
            end = args.pwm_end

            if run == 'input':
                pstart, pend = start, end
            elif run == 'output':
                pstart, pend = utilityFunctions.updateStartEnd(start,
                                                               end, removed)
            else:
                pstart = 0
                pend = len(arr[0, :])

            if end < start:
                print("Error! The start position should be less \
than the end position for the position weight matrix")
                exit()
            log.info("Position matrices will show between positions \
%i and %i" % (start, end))
        else:
            pstart = 0
            pend = len(c_arr[0, :])

        log.info("Generating position frequency matrix")
        if not args.silent:
            print("Generating position frequency matrix")

        subarr = c_arr[:, pstart:pend]
        # Make the position frequency matrix
        PFM, RNA = matrices.makePFM(subarr, typ)
        # Make the second matrix where needed
        if args.pwm_freqtype == 'calc2':
            PFM2 = matrices.makePFM(c_arr, typ)
        else:
            PFM2 = None

        # Calculate the frequency and the alpha matrices
        freq = matrices.getFreq(args.pwm_freqtype, log, typ, RNA, PFM, PFM2)
        alpha = matrices.getAlpha(args.pwm_alphatype, log, PFM, freq,
                                  args.pwm_alphaval)

        log.info("Generating position probability matrix")
        if not args.silent:
            print("Generating position probability matrix")
        # Calculate the PPM from the PFM
        PPM = matrices.makePPM(PFM, alpha=alpha)

        log.info("Generating position weight matrix")
        if not args.silent:
            print("Generating position weight matrix")

        # Calculate the PWM from the PPM
        PWM = matrices.makePWM(PPM, freq)

        # Save all the matrices
        PFM.to_csv("%s_pfm_%s.txt" % (args.outfile_stem, run),  sep="\t")
        PPM.to_csv("%s_ppm_%s.txt" % (args.outfile_stem, run),   sep="\t")
        PWM.to_csv("%s_pwm_%s.txt" % (args.outfile_stem, run),   sep="\t")

        if args.pwm_output_meme:
            # Save the PPM matrix in the format required for MEME input
            # https://meme-suite.org/meme/doc/meme-format.html
            matrices.memeFormat(PPM, typ, RNA, freq, "%s_ppm_meme_%s.txt" % (
                args.outfile_stem, run), args.outfile_stem)

        if args.pwm_output_blamm:
            # Save the PFM matrix in the format required for BLAMM input
            # https://github.com/biointec/blamm
            out = open("%s_pfm_blamm_%s.txt" % (args.outfile_stem, run), "w")
            out.write(">alignment_motif\n")
            for ind, row in zip(PFM.index.values, PFM.values):
                out.write("%s\t[\t%s\t]\n" % (ind,
                                              "\t".join(
                                                  [str(x) for x in row])))
            out.close()
