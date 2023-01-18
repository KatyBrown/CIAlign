#! /usr/bin/env python

import unittest
import pytest
from unittest import mock
from parameterized import parameterized
import logging
import numpy as np
import os
import CIAlign.parsingFunctions as parsingFunctions
import CIAlign.utilityFunctions as utilityFunctions
from tests.helperFunctions import readMSA

class CleaningFunctionsTests(unittest.TestCase):

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/example1.fasta")
        self.relativePositions = list(range(0, len(self.in_array[0])))
        self.rm_file = 'mock_rmfile.txt'
        self.keeps = {'all_rowwise': np.array(['Seq7']),
                      'crop_ends': np.array([]),
                      'remove_divergent': np.array([]),
                      'remove_short': np.array([])}

    def tearDown(self):
        os.remove(self.rm_file)

    @parameterized.expand([
            [0.05, 0.1, "./tests/test_files/crop_ends_cleaned.fasta", ],
            [0.6, 0.1, "./tests/test_files/example1.fasta", ],
            [0.05, 0.05, "./tests/test_files/example1.fasta", ],
    ])
    def testCropEnds(self, mingap, redefine_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, names = parsingFunctions.cropEnds(self.in_array, self.nams, self.relativePositions, self.rm_file, mock_debug, self.keeps, mingap, redefine_perc)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())

    @parameterized.expand([
            [3, 100, 5, 0.5, "./tests/test_files/remove_insertions_cleaned.fasta", ],
            [20, 100, 5, 0.1, "./tests/test_files/example1.fasta", ],
            [3, 100, 30, 0.9, "./tests/test_files/example1.fasta", ],
    ])
    def testRemoveInsertions(self, min_size, max_size, min_flank, 
                             min_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r, self.relativePositions = parsingFunctions.removeInsertions(self.in_array, self.relativePositions, self.rm_file, mock_debug, min_size, max_size, min_flank, min_perc)
        # check if dimensions are equal first
        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())

    @parameterized.expand([
            [0.75, "./tests/test_files/remove_divergent_cleaned.fasta", ],
            [0.1,  "./tests/test_files/example1.fasta", ],
    ])
    def testRemoveDivergent(self, percidentity, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r = parsingFunctions.removeDivergent(self.in_array, self.nams, self.rm_file, mock_debug, self.keeps, percidentity)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(result_ali[:,0].size, exp_array[:,0].size)
        self.assertGreaterEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())

    @parameterized.expand([
            [50, "./tests/test_files/remove_short_cleaned.fasta", ],
            [10,  "./tests/test_files/example1.fasta", ],
    ])
    def testRemoveShort(self, min_length, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r = parsingFunctions.removeTooShort(self.in_array, self.nams, self.rm_file, mock_debug, self.keeps, min_length)

        self.assertTrue((result_ali == exp_array).all())
        self.assertGreaterEqual(len(self.in_array), len(result_ali))

    def testRemoveGapOnly(self):
        expected = "./tests/test_files/remove_gaponly_cleaned.fasta"
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r, self.relativePositions = parsingFunctions.removeGapOnly(self.in_array, self.relativePositions, self.rm_file, mock_debug)

        self.assertTrue((result_ali == exp_array).all())
        self.assertEqual(len(self.in_array), len(result_ali))


    @parameterized.expand([
            [0.5, 0.5, 5, "./tests/test_files/crop_divergent_cleaned.fasta", ],
            [0.15, 0.15, 3,  "./tests/test_files/example1.fasta", ],
            [0.95, 0.15, 3,  "./tests/test_files/crop_divergent_cleaned2.fasta", ],
    ])
    def testCropDivergent(self, min_ident, min_nongap, buffer, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r, self.relativePositions = parsingFunctions.cropDivergent(self.in_array,
                                                                                   self.relativePositions,
                                                                                   self.rm_file,
                                                                                   mock_debug,
                                                                                   min_ident,
                                                                                   min_nongap, 
                                                                                   buffer)
        # check if dimensions are equal first
        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())


class ParsingTestComplex(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/complex_parsing.fasta")
        self.relativePositions = list(range(0, len(self.in_array[0])))
        self.rm_file = 'mock_rmfile.txt'
        self.keeps = {'all_rowwise': np.array(['Seq7']),
                      'crop_ends': np.array([]),
                      'remove_divergent': np.array([]),
                      'remove_short': np.array([])}

    def tearDown(self):
        if os.path.exists(self.rm_file):
            os.remove(self.rm_file)

    @parameterized.expand([
            [1, 100, 5, 0.5, "./tests/test_files/complex_insertions_cleaned.fasta", ],
    ])
    def testRemoveInsertionsComplex(self, min_size, max_size, min_flank, 
                             min_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r, self.relativePositions = parsingFunctions.removeInsertions(self.in_array, self.relativePositions, self.rm_file, mock_debug, min_size, max_size, min_flank, min_perc)
        # check if dimensions are equal first
        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())
    
    @parameterized.expand([
            [0.05, 0.5, "./tests/test_files/complex_crop_ends_cleaned.fasta", ],
    ])
    def testCropEndsComplex(self, mingap, redefine_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, names = parsingFunctions.cropEnds(self.in_array, self.nams, self.relativePositions, self.rm_file, mock_debug, self.keeps, mingap, redefine_perc)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())

    @parameterized.expand([
            [0.01, 0.5, "./tests/test_files/complex_crop_ends_cleaned.fasta", ],
    ])
    def testCropEndsError(self, mingap, redefine_perc, expected):
        exp_array, names = readMSA(expected)
        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, names = parsingFunctions.cropEnds(
                self.in_array, self.nams, self.relativePositions,
                self.rm_file, mock_debug, self.keeps,
                mingap, redefine_perc)
        captured = self.capsys.readouterr()
        self.assertTrue("Given the length of" in captured.out)

    @parameterized.expand([
            [0.75, "./tests/test_files/complex_remove_divergent_cleaned.fasta", ],
    ])
    def testRemoveDivergentComplex(self, percidentity, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r = parsingFunctions.removeDivergent(self.in_array, self.nams, self.rm_file, mock_debug, self.keeps, percidentity)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(result_ali[:,0].size, exp_array[:,0].size)
        self.assertGreaterEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())

class CropEndsTestComplex1(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/uORF_nt.fasta")
        self.relativePositions = list(range(0, len(self.in_array[0])))
        self.rm_file = 'mock_rmfile.txt'
        self.keeps = {'all_rowwise': np.array(['Seq7']),
                      'crop_ends': np.array([]),
                      'remove_divergent': np.array([]),
                      'remove_short': np.array([])}

    def tearDown(self):
        if os.path.exists(self.rm_file):
            os.remove(self.rm_file)

    @parameterized.expand([
            [0.1, 0.001, "./tests/test_files/overlap_crop_ends_cleaned.fasta", ],
    ])
    def testCropEndsOverlap(self, mingap, redefine_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, names = parsingFunctions.cropEnds(self.in_array, self.nams, self.relativePositions, self.rm_file, mock_debug, self.keeps, mingap, redefine_perc)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())


class CropEndsTestComplex2(unittest.TestCase):

    def setUp(self):
        self.in_array, self.nams = readMSA("./tests/test_files/example2.fasta")
        self.relativePositions = list(range(0, len(self.in_array[0])))
        self.rm_file = 'mock_rmfile.txt'
        self.keeps = {'all_rowwise': np.array(['Seq7']),
                      'crop_ends': np.array([]),
                      'remove_divergent': np.array([]),
                      'remove_short': np.array([])}

    def tearDown(self):
        if os.path.exists(self.rm_file):
            os.remove(self.rm_file)

    @parameterized.expand([
            [0.1, 0.01, "./tests/test_files/example2_crop_ends_clean.fasta", ],
    ])
    def testCropEndsNoCrop(self, mingap, redefine_perc, expected):
        exp_array, names = readMSA(expected)

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, names = parsingFunctions.cropEnds(self.in_array, self.nams, self.relativePositions, self.rm_file, mock_debug, self.keeps, mingap, redefine_perc)

        self.assertEqual(result_ali[0,:].size, exp_array[0,:].size)
        self.assertEqual(len(self.in_array), len(result_ali))
        self.assertTrue((result_ali == exp_array).all())


class RemoveTestEmpty(unittest.TestCase):

    def setUp(self):
        self.in_array, self.nams = utilityFunctions.FastaToArray(
            "./tests/test_files/empty.fasta")
        self.relativePositions = list()
        self.rm_file = 'mock_rmfile.txt'
        self.keeps = {'all_rowwise': np.array(['Seq7']),
                      'crop_ends': np.array([]),
                      'remove_divergent': np.array([]),
                      'remove_short': np.array([])}

    def tearDown(self):
        if os.path.exists(self.rm_file):
            os.remove(self.rm_file)

    @parameterized.expand([
            [50],
    ])
    def testRemoveShortEmpty(self, min_length):

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r = parsingFunctions.removeTooShort(self.in_array,
                                                            self.nams,
                                                            self.rm_file,
                                                            mock_debug,
                                                            self.keeps,
                                                            min_length)
        self.assertTrue(np.size(result_ali) == 0)

    def testRemoveGapEmpty(self):
        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            result_ali, r, self.relativePositions = parsingFunctions.removeGapOnly(self.in_array, self.relativePositions, self.rm_file, mock_debug)
        self.assertTrue(np.size(result_ali) == 0)



if __name__ == '__main__':
    unittest.main(warnings='ignore')
