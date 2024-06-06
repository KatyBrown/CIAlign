#! /usr/bin/env python

'''
python3 -m unittest tests.consensusSeqTest
in CIAlign folder
'''

import unittest
from unittest import mock
from parameterized import parameterized
import logging
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import image

import CIAlign
import CIAlign.consensusSeq as consensusSeq
import CIAlign.utilityFunctions as utilityFunctions
from tests.helperFunctions import readMSA
import shutil
import skimage.metrics

class ConsensusSeqTests(unittest.TestCase):

    def testGetAxisUnits(self):
        fig = plt.figure()
        sub = fig.add_subplot(1, 1, 1)
        expected = dict({'axis_height_px': 369,
                         'axis_width_px': 496,
                         'axis_bottom': 0,
                         'axis_top': 1,
                         'axis_left': 0,
                         'axis_right': 1,
                         'axis_width_u': 1,
                         'axis_height_u': 1,
                         'u_height_px': 369,
                         'u_width_px': 496})

        D = consensusSeq.getAxisUnits(sub)
        D_int = dict([key, int(value)] for key, value in D.items())

        self.assertEqual(D_int, expected)

    @parameterized.expand([
            [1, 479, ],
            [0.5, 239, ],
            [7, 3353, ],
    ])
    def testGetFontSize(self, height_u, expected):
        fig = plt.figure()
        sub = fig.add_subplot(1, 1, 1)

        height = consensusSeq.getFontSize(fig, sub, height_u)

        self.assertEqual(int(height), expected)

    @parameterized.expand([
            [75, 300, 18],
            [1000, 300, 240],
            [100, 100, 72],
    ])
    def testPixelsToPoints(self, pixels, dpi, expected):
        points = consensusSeq.PixelsToPoints(pixels, dpi)

        self.assertEqual(int(points), expected)

    @parameterized.expand([
            ["./tests/test_files/consensus_example_nt.fasta", "majority",
             ['G', 'G', 'G', '-', '-', '-', '-', '-', '-', 'U', 'U', 'A', 'U', 'C', 'U', 'C'],
             [1, 0.83, 0.83, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 1, 1, 1, 1, 1, 1, 0.83],],
            ["./tests/test_files/consensus_example_nt.fasta", "majority_nongap",
             ['G', 'G', 'G', 'C', 'U', 'C', 'U', 'U', 'A', 'U', 'U', 'A', 'U', 'C', 'U', 'C'],
             [1, 0.83, 0.83, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 1, 1, 1, 1, 1, 1, 0.83],],
            ["./tests/test_files/consensus_example_aa.fasta", "majority",
            ['L', 'Q', 'N', 'P', 'R', 'V', 'T', 'Q', 'H', 'S', 'V', '-', '-', '-', '-', '-'],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.03, 0.03, 0.03, 0.03, 0.03],],
            ["./tests/test_files/consensus_example_aa.fasta", "majority_nongap",
            ['L', 'Q', 'N', 'P', 'R', 'V', 'T', 'Q', 'H', 'S', 'V', 'P', 'V', 'R', 'R', 'Y'],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.03, 0.03, 0.03, 0.03, 0.03],],
            ["./tests/test_files/consensus_example_aa_gaps.fasta", "majority_nongap",
            ['L', 'Q', 'N', 'P', 'R', 'V', 'T', 'Q', 'H', 'S', 'V', 'P', 'V', 'R', 'N', 'Y'],
            [0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
             0.03, 0.03, 0.03, 0, 0.03],],            
    ])
    def testFindConsensus(self, MSA, type, expected_consensus, expected_coverage):
        alignment, names = readMSA(MSA)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            consensus, coverage = consensusSeq.findConsensus(alignment, logger, type)
        coverage_rounded = [round(num, 2) for num in coverage]
        self.assertEqual(consensus, expected_consensus)
        self.assertEqual(coverage_rounded, expected_coverage)

    @parameterized.expand([
            [{'-': 1, 'C': 1, 'G': 3, 'U': 1}, 6, 'nt',
            {'A': 0.0, 'G': 0.13, 'T': 0.0, 'C': 0.04,
            'N': 0.0, '-': 0.0, 'U': 0.04, 'R': 0.0, 'Y': 0.0, 'S': 0.0,
            'W': 0.0, 'K': 0.0, 'M': 0.0, 'B': 0.0, 'D': 0.0, 'H': 0.0, 'V': 0.0, 'X': 0.0},
            {'A': 0, 'G': 1.56, 'T': 0, 'C': 1.54, 'N': 0, '-': 0, 'U': 1.54,
            'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'X': 0}],

            [{'V': 32}, 32, 'aa',
            {'D': 0.0, 'E': 0.0, 'C': 0.0, 'M': 0.0, 'K': 0.0, 'R': 0.0,
            'S': 0.0, 'T': 0.0, 'F': 0.0, 'Y': 0.0, 'N': 0.0, 'Q': 0.0, 'G': 0.0,
            'L': 0.0, 'V': 3.89, 'I': 0.0, 'A': 0.0, 'W': 0.0, 'H': 0.0, 'P': 0.0,
             'X': 0.0, '-': 0.0, 'B': 0.0, 'Z': 0.0, 'J': 0.0, '*': 0.0, 'U': 0.0, 'O': 0.0},
            {'D': 0, 'E': 0, 'C': 0, 'M': 0, 'K': 0, 'R': 0, 'S': 0, 'T': 0,
             'F': 0, 'Y': 0, 'N': 0, 'Q': 0, 'G': 0, 'L': 0, 'V': 4.32, 'I': 0,
              'A': 0, 'W': 0, 'H': 0, 'P': 0, 'X': 0, '-': 0, 'B': 0, 'Z': 0, 'J': 0, '*': 0, 'U': 0, 'O': 0}],
            
            [{}, 0, 'nt',
             {'A': 0, 'G': 0, 'T': 0, 'C': 0, 'N': 0, '-': 0, 'U': 0,
              'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0,
              'D': 0, 'H': 0, 'V': 0, 'X': 0},
             {'A': 0, 'G': 0, 'T': 0, 'C': 0, 'N': 0, '-': 0, 'U': 0,
              'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0,
              'D': 0, 'H': 0, 'V': 0, 'X': 0}]

    ])
    def test_calc_entropy(self, count, seq_count, type, expected_height, expected_info):
        height, info, freq = consensusSeq.calc_entropy(count, seq_count, type)
        height_rounded = dict([key, round(value, 2)] for key, value in height.items())
        info_rounded = dict([key, round(value, 2)] for key, value in info.items())
        self.assertEqual(height_rounded, expected_height)
        self.assertEqual(info_rounded, expected_info)

class ConsensusSeqNtLettersTests(unittest.TestCase):

    def setUp(self):
        self.coloursNT = CIAlign.utilityFunctions.getNtColours()

    def tearDown(self):
        # for each possible base/aa delete plot
        for base in self.coloursNT.keys():
            b = base.replace("*", "stop")
            os.remove("%s_temp.png" % b)

    def testGetLetters(self):
        consensusSeq.getLetters('nt','monospace',500)

        for base in self.coloursNT.keys():
            b = base.replace("*", "stop")
            self.assertTrue(os.path.isfile("%s_temp.png" % b))

class ConsensusSeqAALettersTests(unittest.TestCase):

    def setUp(self):
        self.coloursAA = CIAlign.utilityFunctions.getAAColours()

    def tearDown(self):
        # for each possible base/aa delete plot
        for base in self.coloursAA.keys():
            b = base.replace("*", "stop")
            os.remove("%s_temp.png" % b)

    def testGetLetters(self):
        consensusSeq.getLetters('aa','monospace',500)

        for base in self.coloursAA.keys():
            b = base.replace("*", "stop")
            self.assertTrue(os.path.isfile("%s_temp.png" % b))

class ConsensusSeqCoveragePlotTest(unittest.TestCase):

    def setUp(self):
        self.dest = 'coverage_test.png'

    def tearDown(self):
        os.remove(self.dest)

    def testMakeLinePlot(self):
        consensusSeq.makeLinePlot([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.03, 0.03, 0.03, 0.03, 0.03], self.dest, "x")
        self.assertTrue(os.path.isfile(self.dest))

class ConsensusSeqPlotResidueFrequenciesTest(unittest.TestCase):

    def setUp(self):
        self.dest = 'resfreq_test.png'

    def tearDown(self):
        os.remove(self.dest)

    @parameterized.expand([['tests/test_files/example1.fasta', 'nt'],
                           ['tests/test_files/example2.fasta', 'aa']])
    def testPlotResidueFrequencies(self, fasta, typ):
        alignment, names = readMSA(fasta)
        expected = image.imread("tests/test_files/expected_resfreq_%s.png" % typ)[500:1000, 500:1000]
        consensusSeq.plotResidueFrequencies(alignment, typ, self.dest)
        imi = image.imread(self.dest)[500:1000, 500:1000]
        simi = skimage.metrics.structural_similarity(expected,
                                                     imi,
                                                     channel_axis=-1,
                                                     data_range=1)
        self.assertTrue(simi > 0.85)


class ConsensusSeqResidueChangeCountTest(unittest.TestCase):

    def setUp(self):
        self.dest = 'changecount_test.png'

    def tearDown(self):
        os.remove(self.dest)

    @parameterized.expand([['tests/test_files/example1.fasta']])
    def testResidueChangeCount(self, fasta):
        alignment, names = readMSA(fasta)
        expected = image.imread("tests/test_files/expected_changecount.png")[500:1000, 500:1000]
        consensusSeq.residueChangeCount(alignment, 'nt', self.dest)
        imi = image.imread(self.dest)[500:1000, 500:1000]
        simi = skimage.metrics.structural_similarity(expected,
                                                     imi,
                                                     channel_axis=-1,
                                                     data_range=1)
        self.assertTrue(simi > 0.85)


class ConsensusSeqSequenceLogoTest(unittest.TestCase):

    def setUp(self):
        self.dest = 'seq_logo_test.png'

    def tearDown(self):
        os.remove(self.dest)

    @parameterized.expand([['tests/test_files/consensus_example_nt.fasta', 'nt', 1, 0, 50],
                           ['tests/test_files/consensus_example_nt.fasta', 'nt', 0, 0, 50],
                           ['tests/test_files/consensus_example_nt.fasta', 'nt', 1, 15, 10],
                           ['tests/test_files/consensus_example_aa.fasta', 'aa', 0, 100, 50],
                           ['tests/test_files/consensus_example_long.fasta', 'aa', 20, 200, 50]])
    def testSequenceLogo(self, fasta, typ, start, end, figrowlength):
        alignment, names = readMSA(fasta)
        consensusSeq.sequence_logo(alignment, self.dest, typ=typ, start=start,
                                   end=end, figrowlength=figrowlength)
        self.assertTrue(os.path.isfile(self.dest))


class ConsensusSeqCoverageSequenceLogoBarTest(unittest.TestCase):

    def setUp(self):
        self.dest = 'seq_logo_bar_test.png'

    def tearDown(self):
        os.remove(self.dest)


    @parameterized.expand([['tests/test_files/consensus_example_nt.fasta', 'nt', 1, 0, 50],
                           ['tests/test_files/consensus_example_nt.fasta', 'nt', 0, 0, 50],
                           ['tests/test_files/consensus_example_nt.fasta', 'nt', 1, 15, 10],
                           ['tests/test_files/consensus_example_aa.fasta', 'aa', 0, 100, 50],
                           ['tests/test_files/consensus_example_long.fasta', 'aa', 20, 200, 50]])
    def testSequenceBarLogo(self, fasta, typ, start, end, figrowlength):
        alignment, names = readMSA(fasta)
        consensusSeq.sequence_bar_logo(alignment, self.dest, typ,
                                       start=start, end=end,
                                       figrowlength=figrowlength)
        self.assertTrue(os.path.isfile(self.dest))

class ConsensusSeqConservation(unittest.TestCase):
    def setUp(self):
        self.arr, nams = utilityFunctions.FastaToArray(
            "tests/test_files/consensus_example_nt.fasta")

    def tearDown(self):
        pass

    @parameterized.expand([[np.array([0.388, 0.224, 0.764, 0.273, 0.273,
                                      0.273, 0.273, 0.273, 0.273, 1.639,
                                      0.989, 0.989, 0.388, 0.388, 0.388,
                                      0.764]),
                            np.array([0.868, 0.95, 0.5, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.451, 0.451, 0.868, 0.868, 0.868,
                                      0.5])]])
    def testCalcConservationAli(self, expected_heights, expected_ents):
        heights, ents = consensusSeq.calcConservationAli(self.arr, 'nt')
        heights = np.array(heights).round(3)
        ents = np.array(ents).round(3)
        self.assertTrue(np.array_equal(heights, expected_heights))
        self.assertTrue(np.array_equal(ents, expected_ents))