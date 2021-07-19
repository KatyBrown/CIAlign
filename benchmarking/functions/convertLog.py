#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
try:
   import CIAlign.utilityFunctions as utilityFunctions
except ImportError:
   import utilityFunctions

def convert(tool, arr, nams, trimfile, logfile, outfile, percs=[]):
    if tool == "trimal":
        convertTrimalLog(arr, nams, trimfile, logfile, outfile)
    elif tool == 'gblocks':
        convertGBlocksLog(arr, nams, trimfile, logfile, outfile)
    elif tool == "zorro":
        convertZorroLog(arr, nams, trimfile, logfile, outfile, percs)


def convertTrimalLog(arr, nams, trimfile, logfile, outfile):
    '''
    Convert the trimal --colnumbering output to resemble the CIAlign
    "removed" file.
    '''
    t_arr, t_nams = utilityFunctions.FastaToArray(trimfile, "")
    width = np.shape(arr)[1]
    t_width = np.shape(t_arr)[1]
    remaining = [int(x.strip())
                 for x in open(logfile).readlines(
                         )[0].strip().split("\t")[1].split(",")]
    remaining = set(remaining)
    all_ints = set(np.arange(0, width))
    removed = sorted(list(all_ints - remaining))
    assert len(removed) == width - t_width
    
    # Check removing these columns gives the trimal output
    new_arr = arr[:, np.array(list(remaining))]
    assert (new_arr == t_arr).all()
    out = open(outfile, "w")
    out.write("other\t%s\n" % (",".join([str(x) for x in removed])))
    out.close()

def convertGBlocksLog(arr, nams, trimfile, logfile, outfile):
    '''
    Convert the GBlocks txt output to resemble the CIAlign
    "removed" file.
    '''
    t_arr, t_nams = utilityFunctions.FastaToArray(trimfile, "")
    width = np.shape(arr)[1]
    t_width = np.shape(t_arr)[1]
    full = ""
    with open(logfile) as infile:
        for line in infile:
            if line.startswith("Gblocks") and not "Results" in line:
                string = line.split(" ")[-1].strip()
                full += string
    removed = set(list(np.where(np.array(list(full)) != "#")[0]))
    kept = np.where(np.array(list(full)) == "#")[0]
    removed = sorted(list(removed))
    assert len(removed) == width - t_width
    
    # Check removing these columns gives the gblocks output
    new_arr = arr[:, kept]
    assert (new_arr == t_arr).all()

    out = open(outfile, "w")
    out.write("other\t%s\n" % (",".join([str(x) for x in removed])))
    out.close()

def convertZorroLog(arr, nams, trimfile, logfile, outfile, thresh):
    '''
    Convert the Zorro output to resemble the CIAlign "removed" file
    '''
    scores = [float(x.strip()) for x in open(logfile).readlines()]
    scores = np.array(scores)
    which = np.where(scores < thresh)[0]
    keeps = np.where(scores >= thresh)[0]
    removed = set(which)
    new_arr = arr[:, keeps]
    utilityFunctions.writeOutfile(trimfile, new_arr, nams, removed)
    out = open(outfile, "w")
    out.write("other\t%s\n" % (",".join([str(x) for x in removed])))
    out.close()
    

def convertGUIDANCELog(arr, nams, trimfile, logfile, outfile):
    '''
    Convert the GUIDANCE output to resemble the CIAlign "removed" file
    '''
    

    trimfile_cols, trimfile_rows, out_trimmed = trimfile
    logfile_cols, logfile_rows = logfile
    
    t_arr_cols, t_nams_cols = utilityFunctions.FastaToArray(trimfile_cols)
    removed_cols = [int(line.strip().split("\t")[0].split(" ")[-1])
                    for line in open(logfile_cols).readlines()]
    removed_cols = np.array(removed_cols) - 1
    all_ints = set(np.arange(0, np.shape(arr)[1]))
    keep = sorted(list(all_ints - set(removed_cols)))
    if os.path.exists("%s.With_Names" % trimfile_rows) and os.path.exists("%s.With_Names" % logfile_rows):
        t_arr_rows, t_nams_rows = utilityFunctions.FastaToArray(
            "%s.With_Names" % trimfile_rows)
        t_arr_rows_rm, t_nams_rows_rm = utilityFunctions.FastaToArray("%s.With_Names" % logfile_rows)        
    elif os.path.exists("%s.With_Names" % trimfile_rows):
        t_arr_rows, t_nams_rows = utilityFunctions.FastaToArray(
            "%s.With_Names" % trimfile_rows)
        t_arr_rows_rm, t_nams_rows_rm = utilityFunctions.FastaToArray(logfile_rows)
    else:
        t_arr_rows, t_nams_rows = np.array([]), list()
        t_arr_rows_rm, t_nams_rows_rm = utilityFunctions.FastaToArray(
            "%s.With_Names" % logfile_rows)
    
    assert len(t_nams_rows) + len(t_nams_rows_rm) == len(nams)
    assert len(removed_cols) + np.shape(t_arr_cols)[1] == np.shape(arr)[1]
    
    assert (arr[:, keep] == t_arr_cols).all()
    
    allnams = sorted(t_nams_rows + t_nams_rows_rm)

    assert allnams == sorted(nams)
    out = open(outfile, "w")
    out.write("other\t%s\n" % (",".join([str(x)
                                         for x in sorted(removed_cols)])))
    out.write("otherc\t%s\n" % (",".join([str(x)
                                          for x in sorted(t_nams_rows_rm)])))
    out.close()
    which_nams = np.where(np.isin(nams, t_nams_rows))[0]

    new_arr = arr[which_nams, ]
    new_arr = new_arr[:, keep]
    utilityFunctions.writeOutfile(out_trimmed, new_arr, nams, t_nams_rows_rm)