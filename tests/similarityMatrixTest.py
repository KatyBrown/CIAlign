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

import pandas as pd

import CIAlign

from CIAlign.similarityMatrix import calculateSimilarityMatrix
from tests.helperFunctions import readMSA

class SimMatrixTests(unittest.TestCase):

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/remove_gaponly_cleaned.fasta")

    def tearDown(self):
        pass

    @parameterized.expand([
        [0, "./tests/test_files/sim_matrix_input_similarity.tsv", ],
        [1, "./tests/test_files/sim_matrix2_input_similarity.tsv", ],
        [2, "./tests/test_files/sim_matrix3_input_similarity.tsv", ],
    ])
    def testSimilarityMatrix(self, keepgaps, expected):
        exp_matrix = pd.read_csv(expected, sep="\t")
        exp_matrix = np.array(exp_matrix)
        exp_matrix = np.delete(exp_matrix, 0, 1)

        minoverlap = 1
        outfile = None
        dp = 4 # does not effect the return value of function but what is saved in the output file

        result_matrix = calculateSimilarityMatrix(self.in_array, self.nams, minoverlap, keepgaps, outfile, dp)

        self.assertTrue((np.round(result_matrix, dp) == exp_matrix).all())
        self.assertEqual(len(self.in_array), len(result_matrix))


if __name__ == '__main__':
    unittest.main(warnings='ignore')
