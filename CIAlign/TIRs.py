#!/usr/bin/env python3
import Bio.Seq as Seq
import numpy as np
import re
try:
    import CIAlign.consensusSeq as consensusSeq
except ImportError:
    import consensusSeq


def find_TIR_pos(consensus, mingap, minlen, maxlen, outstem):
    cons = "".join(consensus)
    consensus_rev = Seq.reverse_complement(cons)
    for L in np.arange(maxlen, minlen, -2):
        for wind in np.arange(0, len(cons)-L, 1):
            substr = cons[int(wind):int(wind+L)]
            matches = re.finditer(substr, consensus_rev)
            i = 0
            for m in matches:
                i += 1
            if i != 0:
                result = next(re.finditer(substr, consensus_rev))
                left_start, left_end = int(wind), int(wind+L)
                matchspan = result.span()
                right_start = len(cons) - matchspan[1]
                right_end = len(cons) - matchspan[0]
                if right_start - left_end > mingap:
                    s1 = cons[left_start:left_end]
                    s2 = cons[right_start:right_end]                 
                    if s1 !=  s2:
                        s2_rev = Seq.reverse_complement(s2)
                        assert s2_rev == s1
                        out = open("%s_aligned_TIRs.fasta" % outstem, "w")
                        out.write(">S1\n%s\n>S2\n%s\n" % (s1, s2_rev))
                        out.close()
                        return (left_start, left_end, right_start, right_end)


def findTIRs(arr, mingap, minlen, maxlen, outstem):
    cons = consensusSeq.findConsensus(arr, "", 'majority_nongap')[0]
    cons = np.array(cons)
    tir_pos = find_TIR_pos(cons, mingap, minlen, maxlen, outstem)
    if tir_pos is None:
        return (None, None, None, None)
    else:
        return (tir_pos)