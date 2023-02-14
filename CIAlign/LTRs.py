#!/usr/bin/env python3
import numpy as np
import Bio.pairwise2 as pairwise2
import itertools
import re
import matplotlib.pyplot as plt
try:
    import CIAlign.consensusSeq as consensusSeq
except ImportError:
    import consensusSeq

def plotMScores(scores, window_int, outfile):
    f = plt.figure(figsize=(5, 5), facecolor='white')
    a = f.add_subplot(111)
    a.imshow(scores, vmin=100, vmax=np.max(scores), cmap='binary')
    a.set_xticks(np.arange(0, np.shape(scores)[0], 50))
    a.set_yticks(np.arange(0, np.shape(scores)[0], 50))
    a.set_xticklabels(np.arange(0, np.shape(scores)[0]*window_int,
                                50*window_int))
    a.set_yticklabels(np.arange(0, np.shape(scores)[0]*window_int,
                                50*window_int))
    f.savefig(outfile, dpi=150, bbox_inches='tight')


def makeWindows(cons, window_size, window_int):
    starts = np.arange(0, len(cons)-window_size, window_int).astype(int)
    ends = starts + window_size
    ints = np.vstack([starts, ends]).astype(int)
    windows = cons[np.apply_along_axis(lambda x: np.arange(*x), 0, ints).T]
    return (windows)


def scorePairs(windows, window_size, window_int):
    scores = np.zeros([len(windows), len(windows)])
    for i, j in itertools.combinations(np.arange(len(windows)), 2):
        seq1 = "".join(windows[i])
        seq2 = "".join(windows[j])
        s1 = i * window_int
        s2 = j * window_int
        if s2 > s1 + window_size:
            score = pairwise2.align.localxx(seq1,
                                            seq2,
                                            one_alignment_only=True,
                                            score_only=True)
            scores[i, j] = score
            scores[j, i] = score
        else:
            scores[i, j] = 0
            scores[j, i] = 0
    return (scores)


def findLTRPos(cons, scores, windows, window_size, window_int, minscore,
               outstem):
    mscore = np.max(scores)

    if mscore >= minscore:
        w = np.where(scores == mscore)
        pos1 = w[0][0]
        pos2 = w[1][0]
    
        seq1 = "".join(windows[pos1])
        seq2 = "".join(windows[pos2])
        ali = pairwise2.align.localxx(seq1, seq2, one_alignment_only=True)
        out = open("%s_aligned_LTRs.fasta" % outstem, "w")
        out.write(">S1\n%s\n>S2\n%s\n" % (ali[0].seqA, ali[0].seqB))
        out.close()
        start = ali[0][-2]
        end = ali[0][-1]
        s1 = ali[0][0][start:end].replace("-", "")
        s2 = ali[0][1][start:end].replace("-", "")

        window1_start = pos1 * window_int
        window1_end = (pos1 * window_int) + window_size

        window2_start = pos2 * window_int
        window2_end = (pos2 * window_int) + window_size

        orig_s1 = "".join(cons[int(window1_start):int(window1_end)])
        orig_s2 = "".join(cons[int(window2_start):int(window2_end)])

        span1 = re.search(s1, orig_s1).span()
        span2 = re.search(s2, orig_s2).span()


        window1_start_final = window1_start + span1[0]
        window1_end_final = window1_start + (span1[1] - span1[0])

        window2_start_final = window2_start + span2[0]
        window2_end_final = window2_start + (span2[1] - span2[0])
        
        return (window1_start_final, window1_end_final,
                window2_start_final, window2_end_final)


def findLTRs(arr, outstem, window_size, window_int, minscore):
    cons = consensusSeq.findConsensus(arr, "", 'majority_nongap')[0]
    cons = np.array(cons)
    windows = makeWindows(cons, window_size, window_int)
    scores = scorePairs(windows, window_size, window_int)
    outplot = "%s_LTRplot.png" % outstem
    plotMScores(scores, window_int, outplot)
    ltr_pos = findLTRPos(cons, scores, windows, window_size, window_int,
                         minscore, outstem)
    if ltr_pos is None:
        return (None, None, None, None)
    else:
        return (ltr_pos)