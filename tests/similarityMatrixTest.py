#! /usr/bin/env python

import unittest
from parameterized import parameterized
import numpy as np
import os
import pandas as pd
from CIAlign.similarityMatrix import calculateSimilarityMatrix
from tests.helperFunctions import readMSA

class SimMatrixTests(unittest.TestCase):

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/remove_gaponly_cleaned.fasta")
        self.outfile = "tests/test_files/testing_simmatrix.txt"

    def tearDown(self):
        os.remove(self.outfile)

    @parameterized.expand([
        [0, 1, "./tests/test_files/sim_matrix_input_similarity.tsv"],
        [1, 1, "./tests/test_files/sim_matrix2_input_similarity.tsv"],
        [2, 1, "./tests/test_files/sim_matrix3_input_similarity.tsv"],
        [2, 100, "./tests/test_files/sim_matrix4_input_similarity.tsv"]
        ])
    def testSimilarityMatrix(self, keepgaps, minoverlap, expected):
        dp = 4 # does not effect the return value of function but what is saved
        # in the output file

        exp_matrix = pd.read_csv(expected, sep="\t").round(dp)
        exp_matrix = np.array(exp_matrix)
        exp_matrix = np.delete(exp_matrix, 0, 1)
        result_matrix = calculateSimilarityMatrix(self.in_array, self.nams,
                                                  minoverlap, keepgaps,
                                                  dp=dp, outfile=self.outfile)
        
        self.assertTrue((np.round(result_matrix, dp) == exp_matrix).all())
        self.assertEqual(len(self.in_array), len(result_matrix))
        ff = pd.read_csv(self.outfile, sep="\t").round(dp)
        ff = np.array(ff)
        ff = np.delete(ff, 0, 1)
        self.assertTrue(np.array_equal(ff, exp_matrix))


if __name__ == '__main__':
    unittest.main(warnings='ignore')
