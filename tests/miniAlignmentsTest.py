#! /usr/bin/env python

'''
python3 -m unittest tests.miniAlignmentsTest
in CIAlign folder
'''

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
import matplotlib.pyplot as plt

from matplotlib import image

import CIAlign
import CIAlign.miniAlignments as miniAlignments
from tests.helperFunctions import readMSA

class MiniAlignmentsTests(unittest.TestCase):

    def testArrNumeric(self):
        alignment, names = readMSA("./tests/test_files/consensus_example_nt.fasta")
        arr_expected = [[1, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 0, 4, 2, 4, 3],
                    [1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 4, 0, 4, 2, 4, 2],
                    [1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 4, 0, 0, 1, 2, 4],
                    [1, 1, 1, 2, 4, 2, 4, 4, 0, 4, 4, 0, 4, 2, 4, 2],
                    [4, 2, 1, 3, 3, 3, 3, 3, 3, 4, 4, 0, 4, 2, 4, 2],
                    [2, 4, 0, 3, 3, 3, 3, 3, 3, 4, 2, 1, 2, 4, 0, 2],]

        arr_int, colour_map = miniAlignments.arrNumeric(alignment,'nt')

        self.assertTrue((arr_int == arr_expected).all())

    def testDrawMarkUp(self):

        self.assertEqual(0, 0)

    def testDrawMarkUpLegend(self):

        self.assertEqual(0, 0)

class MiniAlignmentsDrawTest(unittest.TestCase):

    def setUp(self):
        self.dest = "./tests/test_files/test_mini.png"
        self.alignment, self.names = readMSA("./tests/test_files/example1.fasta")

    def tearDown(self):
        os.remove("./tests/test_files/test_mini_legend.png")
        os.remove(self.dest)

    def testDrawMiniAlignment(self):
        expected = image.imread('./tests/test_files/expected_mini_ali.png')
        markup_dict = {'remove_divergent': {'Seq1'},
                        'remove_gaponly': {89, 90, 91, 92, 93, 94, 95},
                        'remove_insertions': {22, 23, 24, 25, 26, 27},
                        'crop_ends': {'Seq5': ((np.array([]), np.array([92, 93, 94, 95])))},
                        'remove_short': {'Seq6'}}

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            miniAlignments.drawMiniAlignment(self.alignment, self.names, logger, self.dest,
                                                                'nt', 300, None, 5, 3,
                                                                True, markup_dict, False)

        mini_alignment = image.imread(self.dest)

        self.assertTrue((mini_alignment == expected).all())


# https://www.pluralsight.com/guides/importing-image-data-into-numpy-arrays
# i think PIL can load images from a flat file and convert to np array
# and then numpy.array_equal to check they are the same
# we already depend on PIL because matplotlib needs it
