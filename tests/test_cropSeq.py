#! /usr/bin/env python

import unittest
from unittest import mock
from mock import patch
from parameterized import parameterized, parameterized_class

import sys
import logging
import numpy as np
from Bio import AlignIO
import os
from os import path

import CIAlign

import CIAlign.cropSeq as cropSeq

class cropSeqsTests(unittest.TestCase):

    # todo: more test cases
    @parameterized.expand([
            [0.05, 0.1, '--UC----UCUCUCUCGCGUGUGUGAAAAAA----AAAUUUU------------A', 8, 42],
            [0.05, 0.1, '--UC--AA-----UCUCUCUCGCGUGUGUGAAAAAA----AAAUUUU------------A', 6, 47],
            [0.05, 0.2, '--UC--AA-----UCUCUCUCGCGUGUGUGAAAAAA----AAAUUUU------------A', 13, 47],
    ])
    def test_determineStartEnd(self, mingap, redefine_perc, input, expected_start, expected_end):

        seqs = []
        seqs.append([s for s in input])
        input = np.array(seqs[0])
        print(input)
        start, end = cropSeq.determineStartEnd(input, mingap, redefine_perc)
        print(start, end)

        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)

# todo: test other functions
