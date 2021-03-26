#! /usr/bin/env python

import numpy as np
from Bio import AlignIO

# helper function to read MSA from file into np array
def readMSA(input):
    in_array = []
    names = []
    input_handle = open(input, 'r')
    input_ali = AlignIO.read(input_handle, "fasta")
    for record in input_ali:
        in_array.append(record.seq)
        names.append(record.id)
    in_array = np.array(in_array)
    input_handle.close()

    return in_array, names
