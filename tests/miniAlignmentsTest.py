#! /usr/bin/env python

'''
python3 -m unittest tests.miniAlignmentsTest
in CIAlign folder
'''
from parameterized import parameterized
import unittest
from unittest import mock
import logging
import numpy as np
import os
from matplotlib import image
import CIAlign.miniAlignments as miniAlignments
from tests.helperFunctions import readMSA
import skimage.metrics
import shutil


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

class DrawMarkUpLegendTest(unittest.TestCase):

    def setUp(self):
        self.dest = './tests/test_files/legend'
        self.result = './tests/test_files/legend_legend.png'

    def tearDown(self):
        os.remove(self.result)

    def testDrawMarkUpLegend(self):
        miniAlignments.drawMarkUpLegend(self.dest)

        self.assertTrue(os.path.isfile(self.result))

class MiniAlignmentsDrawTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        if os.path.isfile(self.legend):
            os.remove(self.legend)
        os.remove(self.dest)

    def testDrawMiniAlignmentMarkUp(self):
        self.alignment, self.names = readMSA("./tests/test_files/example1.fasta")
        self.dest = "./tests/test_files/test_mini_markup.png"
        self.legend = "./tests/test_files/test_mini_markup_legend.png"
        expected = image.imread('./tests/test_files/expected_mini_ali_markup.png').round(3)
        markup_dict = {'remove_divergent': {'Seq1'},
                       'remove_gap_only': {89, 90, 91, 92, 93, 94, 95},
                       'remove_insertions': {22, 23, 24, 25, 26, 27},
                       'crop_ends': {'Seq5': ((np.array([1, 2, 3]), np.array([92, 93, 94, 95])))},
                       'remove_short': {'Seq6'},
                       'crop_divergent': {0, 1, 88, 89},
                       'user': {19, 20}}

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            miniAlignments.drawMiniAlignment(self.alignment, self.names,
                                             logger, self.dest,
                                             'nt', 'standard', 300, None, 5, 3,
                                             True, markup_dict, False)

        mini_alignment = image.imread(self.dest).round(3)

        # added a bit of leeway to allow for images created on different
        # machines - they are visually identical but have minor differences
        # in rendering - look for 95% structural similarity
        simi = skimage.metrics.structural_similarity(expected,
                                                     mini_alignment,
                                                     channel_axis=-1,
                                                     data_range=1)

        self.assertTrue(simi > 0.95)

    @parameterized.expand([
            ['./tests/test_files/example1.fasta', './tests/test_files/expected_mini_ali.png', 'nt', True, 'standard', None, 'CBS', None],
            ['./tests/test_files/example2.fasta', './tests/test_files/expected_mini_ali_aa.png', 'aa', False, 'standard', None, 'Bright', None],
            ['./tests/test_files/example1.fasta', './tests/test_files/expected_identity_nt.png', 'nt', True, 'identity', None, 'bone', None],
            ['./tests/test_files/example2.fasta', './tests/test_files/expected_identity_aa.png', 'aa', True, 'identity', None, 'terrain', None],
            ['./tests/test_files/example1.fasta', './tests/test_files/expected_similarity_nt.png', 'nt', True, 'similarity', 'NUC.4.4', 'bone', 'white'],
            ['./tests/test_files/example2.fasta', './tests/test_files/expected_similarity_aa.png', 'aa', True, 'similarity', 'BLOSUM62', 'bone', 'white'],
            ['./tests/test_files/example1.fasta', './tests/test_files/expected_similarity_nt_NUC.4.2.png', 'nt', True, 'similarity', 'NUC.4.2', 'PuOr', 'red'],
            ['./tests/test_files/example2.fasta', './tests/test_files/expected_similarity_aa_PAM10.png', 'aa', True, 'similarity', 'PAM10', 'rainbow', 'blue'],
    ])
    def testDrawMiniAlignment(self, input, expected, type, keep_numbers, plot_type, matrix, pal, gapcol):
        self.alignment, self.names = readMSA(input)
        self.dest = "./tests/test_files/test_mini.png"
        self.legend = ""
        expected = image.imread(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            miniAlignments.drawMiniAlignment(self.alignment, self.names, logger, self.dest,
                                                                type, plot_type, 300,
                                                                title="test",
                                                                width=5,
                                                                height=3,
                                                                markup=False,
                                                                markupdict=None,
                                                                ret=False,
                                                                orig_nams=['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5', 'Seq6', 'Seq7'],
                                                                keep_numbers=keep_numbers,
                                                                force_numbers=False,
                                                                palette=pal,
                                                                plot_identity_palette=pal,
                                                                plot_similarity_palette=pal,
                                                                plot_identity_gap_col=gapcol,
                                                                plot_similarity_gap_col=gapcol,
                                                                sub_matrix_name=matrix)

        mini_alignment = image.imread(self.dest)
        # added a bit of leeway to allow for images created on different
        # machines - they are visually identical but have minor differences
        # in rendering - look for 95% structural similarity
        simi = skimage.metrics.structural_similarity(expected,
                                                     mini_alignment,
                                                     channel_axis=-1,
                                                     data_range=1)
        self.assertTrue(simi > 0.85)


class DrawMarkUpTest(unittest.TestCase):

    def setUp(self):
        self.dest = "./tests/test_files/test_mini.png"
        self.alignment, self.names = readMSA("./tests/test_files/example1.fasta")
        self.ali_width = len(self.alignment[0])
        self.ali_height = len(self.alignment)
        self.markup_dict = {'remove_divergent': {'Seq1'},
                        'remove_gap_only': {89, 90, 91, 92, 93, 94, 95},
                        'remove_insertions': {22, 23, 24, 25, 26, 27},
                        'crop_ends': {'Seq5': ((np.array([]), np.array([92, 93, 94, 95])))},
                        'remove_short': {'Seq6'}}
        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            # make mini plot without markup
            self.mini_plot = miniAlignments.drawMiniAlignment(self.alignment, self.names, logger, self.dest,
                                                                'nt', 'standard', 300, None, 5, 3,
                                                                False, None, True)
            # make mini plot with markup
            self.mini_plot_expected = miniAlignments.drawMiniAlignment(self.alignment, self.names, logger, self.dest,
                                                                'nt', 'standard', 300, None, 5, 3,
                                                                True, self.markup_dict, True)

    def tearDown(self):
        if os.path.isfile("./tests/test_files/test_mini_legend.png"):
            os.remove("./tests/test_files/test_mini_legend.png")
        os.remove(self.dest)

    def testDrawMarkUp(self):
        # add markup to mini plot
        miniAlignments.drawMarkUp(self.mini_plot.axes[0], self.markup_dict, self.names, self.ali_width, self.ali_height)
        # convert plots to numpy arrays to compare
        self.mini_plot.canvas.draw()
        self.mini_plot_expected.canvas.draw()
        data = np.frombuffer(self.mini_plot.canvas.tostring_rgb(), dtype=np.uint8)
        data2 = np.frombuffer(self.mini_plot_expected.canvas.tostring_rgb(), dtype=np.uint8)
        self.assertTrue((data == data2).all())
