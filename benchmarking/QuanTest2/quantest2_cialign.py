#!/usr/bin/env python3
__author__ = "fabian"
__date__ = "2019"

# Adapted from QuanTest2 available via Higgins et al. (2019)
# DOI: 10.1093/bioinformatics/btz552,
# licensed under Creative Commons Attribution License (http://creativecommons.org/licenses/by/4.0/)

import numpy as np
import copy
import itertools
import sys
import ete3
import subprocess

import numpy as np
from Bio import AlignIO

GAPCHARS = ['-', '.']
GAPTHRESH = 0.90
REF_TARGET = 3 # exepect every reference file to have that many references
verbose =  True  # minimum verbosity
VERBOSE =  True #True ## default should be 'False'
VVERBOSE = False #True ## default *definitely* should be 'False'
RUTHLESS = False #False #True ## tidy up files down-loaded from Jpred, should be 'True'
UNKNOWN = -1

nFiles_3 = 0
alnFiles = []
ssFiles  = []
ciaFiles = []
emails = []
cnt = 0
indiScores = []
predicted = []

myBasename=sys.argv[0].rsplit('/',1)[0]
if VERBOSE:
    print("myBasename = {}".format(myBasename))

LOGFILE = "quantest2.log"
try:
    logFile = open(LOGFILE,"w")
except:
    print("OutputError: Could not open file {} for writing".format(LOGFILE))


# helper function to read MSA from file into np array
def readMSA(infile, log=None, outfile_stem=None):
    formatErrorMessage = "The MSA file needs to be in FASTA format."
    nams = []
    seqs = []
    nam = ""
    seq = ""
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue  # todo: test!
            if line[0] == ">":
                seqs.append([s.upper() for s in seq])
                nams.append(nam.upper())
                seq = []
                nam = line.replace(">", "")
            else:
                if len(nams) == 0:
                    if log:
                        log.error(formatErrorMessage)
                    print(formatErrorMessage)
                    exit()
                seq += list(line.replace('-', ''))
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam.upper())
    arr = np.array(seqs[1:])
    return (arr, nams[1:])

#### Commandline() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# put commandline arguments into proper arrays/variables
# return nFiles_3, number of alignment files (= number of secondary structure files)
def Commandline(argv):

    usage = "Usage: "+argv[0]+" <email> <alignment-file1> ... <secondary-file1> ..."

    argaux = []


    for a in argv[1:]:
        if (a.find('@') != -1):
            emails.append(a)
        else:
            argaux.append(a)

    if VVERBOSE:
        print(argaux)

    if len(emails) < 1:
        print("no email address supplied")
        print(usage)
        quit()
    nFiles = len(argaux)
    nFiles_3 = int(nFiles / 3)
    if nFiles < 3:
        print("Commadline Error: at least 2 input files expected")
        print(usage)
        quit()
    if nFiles % 3 != 0:
        print("equal number of alignment files and CIAlign files and secondary structure files expected")
        print(usage)
        quit()
    # 1st half of input files are alignments, 2nd half of input files should be secondary structures
    for i in range(nFiles_3):
        alnFiles.append(argaux[i])
        ssFiles.append(argaux[i+nFiles_3])
        ciaFiles.append(argaux[i+2*nFiles_3])

    print(alnFiles)
    print(ssFiles)
    print(ciaFiles)
    return nFiles_3
# end of def Comanndline() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#### ReadSS() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# read file with secondary structure information
# return NOOFREF, number of secondary structures
def ReadSS(ssFilePtr, refLabel, refSeq):

    NOOFREF = UNKNOWN
    for line_ in ssFilePtr:
        line = line_.rstrip()
        if line == "":
            continue
        if line[0] == ">":
            NOOFREF = NOOFREF+1
            refLabel.append(line[1:])
            refSeq.append("")
        else:
            refSeq[NOOFREF] =  refSeq[NOOFREF]+line
    ssFilePtr.close()
    NOOFREF = NOOFREF+1

    if VVERBOSE:
        print(refSeq)

    return NOOFREF

# end of ReadSS() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#### def ReadAln() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# read alignment file
# return number of sequences in alignment file
def ReadAln(alnFilePtr, label, seq, cnt, NOOFREF, refLabel, REFPOS, family):

    cnt = UNKNOWN
    for line_ in alnFilePtr:
        line = line_.rstrip()
        if line == "":
            continue
        if line[0] == ">":
            cnt = cnt+1
            label.append(line[1:])
            seq.append("")
            for i in range(NOOFREF):
                if label[cnt] == refLabel[i]:
                    REFPOS[i] = cnt
        else:
            if cnt == UNKNOWN:
                print("Format Error: no valid sequence label found")
                quit()
            seq[cnt] = seq[cnt]+line

    cnt = cnt+1

    if VVERBOSE:
        for i in range(cnt):
            print("seq[{}] = {}\t{}".format(i, label[i],seq[i]))


    # check that all reference sequences were present
    for i in range(NOOFREF):
        if REFPOS[i] == UNKNOWN:
            print("WARNING: reference sequence {} ({}) not in alignment".format(i, refLabel[i]))

    # check all sequence lines have the same length
    alnLength = len(seq[0])
    if VVERBOSE:
        print("1st sequence is {} long, there are {} sequences".format(alnLength,cnt))
    for i in range(cnt):
        if len(seq[i]) != alnLength:
            print("WARNING: sequence {} ({}) has length {}, not same as previous sequence/s {}".format(i,label[i],len(seq[i]),alnLength))
            quit()
    if VVERBOSE:
        logFile.write("all aligned sequences in {} have correct lengths\n".format(family))

    return cnt

# end of ReadAln() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



#### def Gblocks() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# put sequence to be predicted on top, remove columns in alignment where 1st sequences has gaps
# write this construct to file, return name of this file
def Gblocks(family, ri, refLabel, refSeq, REFPOS, cnt, seq):

    alnLength = len(seq[0])

    blockfile=family+"-"+str(ri)+".blk"
    try:
        blk = open(blockfile, 'w')
    except:
        print("Could not open file {} for writing (permission?)".format(blockfile))
        quit()

    # print sequence that should be predicted
    if VVERBOSE:
        print("ref {}: {}, position {}, write to {}".format(ri,refLabel[ri],REFPOS[ri],blockfile))
    printline = ""
    for ci in range(alnLength):
        c = seq[REFPOS[ri]][ci]
        # if c != '-' and c != '.':
        printline += c
    # check that the number of residues equals number of ss states
    if len(printline) != len(refSeq[ri]):
        print("number of residues ({}) for reference sequence {} ({}[{}]) does not equal #SS states ({})".format(len(printline),refLabel[ri],family,ri,len(refSeq[ri])))
        # quit()
    elif VVERBOSE:
        logFile.write("reference sequence {}/{} has correct number of residues/states\n".format(refLabel[ri],family))
    blk.write(">{}\n{}\n".format(refLabel[ri],printline))
    if VVERBOSE:
        print(printline)
    # do the other sequences
    k=1
    for j in range(cnt):
        if j == REFPOS[ri]:
            continue
        printline = ""
        for ci in range(alnLength):
            cr = seq[REFPOS[ri]][ci]
            # if cr != '-' and cr != '.':
            printline += seq[j][ci]
        gapCnt = 0
        for gc in GAPCHARS:
            gapCnt = gapCnt + printline.count(gc)
        if (gapCnt / len(refSeq[ri])) > GAPTHRESH:
            continue
        blk.write(">{}\n{}\n".format(k,printline))
        k = k+1
        if VVERBOSE:
            print("{} ({})".format(printline,j))
    blk.close()

    return blockfile

# end of Gblocks() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# extract and score predictions
#### def ExtractAndScore() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# extract jnet file from prediction, retrieved from jpred
def ExtractAndScore(blkReg, jobReg, nmeReg, refReg, indiScores, predicted, cia_check, logFile):

    for j in range(len(jobReg)):

        if VERBOSE:
            print(jobReg[j])
        jnet = jobReg[j]+".jnet"
        tarLog = subprocess.getoutput(["tar -xvf ./"+jobReg[j]+"/"+jobReg[j]+".tar.gz "+jnet])
        #print(tarLog)

        try:
            jnetPtr = open(jnet, "r")
        except FileNotFoundError:
            print("Problem reading file {} ({}, downloaded/un-tar-ed) for reading".format(jnet, blkReg[j]))
            quit()

        for line_ in jnetPtr:
            line = line_.rstrip()
            #print(line)
            colonSplit = line.split(':')
            if colonSplit[0] == "jnetpred":
                states = colonSplit[1].split(',')
            else:
                continue

        if VVERBOSE:
            print("{} has {} states (job {}, seq {})".format(states,len(states),jobReg[j],nmeReg[j]))
            print("{} has {} states".format(refReg[j],len(refReg[j])))

        jnetPtr.close()
        if RUTHLESS:
            subprocess.getoutput(["rm "+jnet])
            subprocess.getoutput(["rm -r "+jobReg[j]])

        aux = "".join(states)
        predicted[j] = aux


        match = 0
        cia_adjust = 0
        ####
        for c in range(len(refReg[j])):

            if(cia_check[j][c] == '!'):
                cia_adjust = cia_adjust + 1
            else:
                r = refReg[j][c] # refReg is the list of all reference structures
                s = states[c - cia_adjust] # result structure from jpred

                if   (r == 'H' or r == 'G' or r == 'I') and (s == 'H' or s == 'G' or s == 'I'):
                    match = match + 1
                elif (r == 'E' or r == 'B')             and (s == 'E' or s == 'B'):
                    match = match + 1
                elif (r == '-' or r == 'C')             and (s == '-' or s == 'C'):
                    match = match + 1
                elif (r != 'H' and r != 'G' and r != 'I' and r != 'E' and r != 'B' and r != '-' and r != 'C') or (s != 'H' and s != 'G' and s != 'I' and s != 'E' and s != 'B' and s != '-' and s != 'C'):
                    print("WARNING: unknown SS status, c={}: r = {} / s = {}".format(c,r,s))

        perc =  match/(len(refReg[j])-cia_adjust)*100
        indiScores[j] = perc
        if VERBOSE:
            print("there are {} matches out of {} = {}%".format(match,len(refReg[j])-cia_adjust,perc))

        logFile.write("{}\t{}\t{}\n".format(nmeReg[j], perc, aux))

    predicted[:] = list(predicted)
    indiScores[:] = list(indiScores)


# end of ExtractAndScore() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



def main():

    nFiles_3 = Commandline(sys.argv)

    if VERBOSE:
        print(emails)
        print(alnFiles)
        print(ssFiles)
    if verbose:
        print("parsed commandline")

    # go through different alignment/secondary files
    eCnt=0
    blkReg = []
    jobReg = []
    nmeReg = []
    refReg = []

    for II in range(nFiles_3):

        alignment = alnFiles[II]
        aux = alignment[0:alignment.rfind(".")]
        family = aux[aux.rfind("/")+1:]
        ylimaf = family.replace('-','_')
        try:
            alnFilePtr = open(alignment,"r")
        except FileNotFoundError:
            print("Input Error: alignment file {} (index {}) does not exist".format(alignment,II))
            quit()


        secondary = ssFiles[II]
        cia_file = ciaFiles[II]
        cia_read, cia_names = msa_replaced, names_replaced = readMSA(cia_file)

        try:
            ssFilePtr = open(secondary,"r")
        except FileNotFoundError:
            print("Input Error: alignment file {} (index {}) does not exist".format(secondary,II))
            quit()

        # refLabel/refSeq are different for every file, initialise inside the main loop
        refLabel = []
        refSeq = []
        NOOFREF = ReadSS(ssFilePtr, refLabel, refSeq)
        if NOOFREF != REF_TARGET:
            print("WARNING: SS-file {} has {} references, expecting {} -- possible problem".format(II, NOOFREF, REF_TARGET))

        refReg.extend(refSeq) # refReg is the list of all reference structures that is kept till the end

        REFPOS = []
        cia_check = []
        for i in range(NOOFREF):
            REFPOS.append(UNKNOWN)

        print(refLabel)
        for lalabel in refLabel:
            for namnam in range(len(cia_names)):
                if cia_names[namnam] == lalabel.upper():
                    cia_check.append(cia_read[namnam])

        label = []
        seq = []
        alnLen = UNKNOWN
        # read fasta files, store as vienna
        cnt = UNKNOWN
        cnt = ReadAln(alnFilePtr, label, seq, cnt, NOOFREF, refLabel, REFPOS, family)

        if VVERBOSE:
            print("refLabel = {}".format(refLabel))
            print("refSeq = {}".format(refSeq))
            print("#r = {}".format(NOOFREF))
            print("refpos = {}".format(REFPOS))
            print("r1 = {}".format(seq[REFPOS[0]]))
            print("r2 = {}".format(seq[REFPOS[1]]))
            print("r3 = {}".format(seq[REFPOS[2]]))

        for ri in range(NOOFREF):
            # do gblocks
            blockfile = Gblocks(family, ri, refLabel, refSeq, REFPOS, cnt, seq)
            myname=ylimaf+"__"+str(ri)
            useE = eCnt % len(emails)
            # this is the most important section of the code -- we actually submit a (modified) alignment to Jpred
            job = ""
            while not job or job.isspace():
                print(blockfile)
                job = subprocess.getoutput(["perl jpredapi submit mode=msa format=fasta email="+emails[useE]+" file="+blockfile+" name="+myname+" silent"])
            if VERBOSE:
                print(job)
            jobList = job.split()
            if jobList[0] == "ERROR:":
                continue
            thisJob = jobList[-1]
            jobReg.append(thisJob)

            if RUTHLESS:
                subprocess.getoutput(["rm "+blockfile])

            eCnt = eCnt + 1
            blkReg.append(blockfile)
            nmeReg.append(myname)

        if VVERBOSE:
            print("blk = {}".format(blkReg))
            print("nme = {}".format(nmeReg))

    if verbose:
        print("submitted gblocks")

    # this is the second most important section of the code
    for j in jobReg:
        retrievalLog = subprocess.getoutput(["perl jpredapi status jobid="+j+" getResults=yes checkEvery=30 silent"])
        if VERBOSE:
            print(retrievalLog)

    if verbose:
        print("retrieved predictions")

    indiScores = [0]*len(jobReg)
    predicted  = [""]*len(jobReg)
    ExtractAndScore(blkReg, jobReg, nmeReg, refReg, indiScores, predicted, cia_check, logFile)

    if verbose:
        print("extracted scores")

    for j in range(len(jobReg)):
        if j % NOOFREF == 0:
            mysum = 0.00
        mysum = mysum + indiScores[j]
        if j % NOOFREF == NOOFREF-1:
            print("{}\t{}".format(alnFiles[int(j/NOOFREF)],mysum/NOOFREF))
            logFile.write("{}\t{}\n".format(alnFiles[int(j/NOOFREF)],mysum/NOOFREF))

if __name__ == "__main__":
    main()
