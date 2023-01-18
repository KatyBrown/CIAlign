#! /usr/bin/env python
import logging
import unittest
import pytest
from unittest import mock
from parameterized import parameterized
import sys
import numpy as np
import os
import warnings
import CIAlign.utilityFunctions as utilityFunctions
from tests.helperFunctions import readMSA


class UtilityFunctionsgetColoursTests(unittest.TestCase):

    def testGetAAColours(self):
        AAcolours = utilityFunctions.getAAColours()
        self.assertEqual(len(AAcolours), 28)
        AAcolours = utilityFunctions.getAAColours(palette='bright')
        self.assertEqual(len(AAcolours), 28)

    def testGetNtColours(self):
        Ntcolours = utilityFunctions.getNtColours()
        self.assertEqual(len(Ntcolours), 18)
        Ntcolours = utilityFunctions.getNtColours(palette='bright')
        self.assertEqual(len(Ntcolours), 18)

    def testGetMarkupColours(self):
        Mcolours = utilityFunctions.getMarkupColours()
        self.assertEqual(len(Mcolours), 7)


class UtilityFunctionsMSAInputTests(unittest.TestCase):

    def setUp(self):
        self.input = "./tests/test_files/example1.fasta"
        self.in_array, self.nams = readMSA(self.input)
        self.rmfile = 'mock_rmfile.txt'

    def tearDown(self):
        if os.path.exists(self.rmfile):
            os.remove(self.rmfile)

    def testReplaceUbyT(self):
        result_ali = utilityFunctions.replaceUbyT(self.in_array, rev=True)
        findU = np.where(result_ali == "U")
        findu = np.where(result_ali == "u")

        self.assertFalse(len(findU[0]) > 0)
        self.assertFalse(len(findU[1]) > 0)
        self.assertFalse(len(findu[0]) > 0)
        self.assertFalse(len(findu[1]) > 0)

    def testReplaceTbyU(self):
        result_ali = utilityFunctions.replaceUbyT(self.in_array, rev=False)
        findT = np.where(result_ali == "T")
        findt = np.where(result_ali == "t")

        self.assertFalse(len(findT[0]) > 0)
        self.assertFalse(len(findT[1]) > 0)
        self.assertFalse(len(findt[0]) > 0)
        self.assertFalse(len(findt[1]) > 0)

    def testUnAlign(self):
        result_ali = utilityFunctions.unAlign(self.in_array)
        findGap = np.where(result_ali == "-")

        self.assertFalse(len(findGap[0]) > 0)
        self.assertFalse(len(findGap[1]) > 0)

    def testFastaToArray(self):
        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            ali, nams = utilityFunctions.FastaToArray(self.input, logger)

        # self.assertEqual(nams.size, self.nams.size)
        self.assertEqual(ali[0,:].size, self.in_array[0,:].size)
        self.assertEqual(len(self.in_array), len(ali))
        self.assertEqual(len(nams), len(self.nams))
        self.assertTrue((ali == self.in_array).all())
        self.assertTrue(nams == self.nams)

    @parameterized.expand([[[22, 23, 24, 25, 26, 27],
                            list(np.arange(0, 98))]])
    def testRemoveColumns(self, rmAbsolute, relativePositions):
        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            arr, nams = utilityFunctions.FastaToArray(self.input, logger)
        newarr, newrel, newrm = utilityFunctions.removeColumns(
            rmAbsolute, relativePositions, arr,
            logger, self.rmfile,'func')
        expected_arr, expected_nams = utilityFunctions.FastaToArray(
            "tests/test_files/remove_insertions_cleaned.fasta")
        self.assertTrue((newarr == expected_arr).all())
        self.assertEqual(sorted(list(newrm)), [22, 23, 24, 25, 26, 27])


class UtilityFunctionsMSAEmpty(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys
        
    def setUp(self):
        self.input = "./tests/test_files/empty.fasta"

    def tearDown(self):
        pass

    @parameterized.expand([
            ['./tests/test_files/corrupt.fasta']])
    def testFastaToArrayCorrupt(self, ali_path):
        logger = logging.getLogger('path.to.module.under.test')

        with self.assertRaises(SystemExit) as cm:
            ali, nams = utilityFunctions.FastaToArray(ali_path, logger)
        captured = self.capsys.readouterr()

        self.assertEqual(cm.exception.code, 1)
        self.assertTrue('needs to be in FASTA format' in captured.out)


class UltilityFunctionsWriteOutfileTest(unittest.TestCase):

    def setUp(self):
        self.input = "./tests/test_files/example1.fasta"
        self.in_array, self.nams = readMSA(self.input)
        self.removed = set()
        self.outfile = "writeOutfile_test.txt"


    def tearDown(self):
        os.remove(self.outfile)

    def testWriteOutfile(self):
        utilityFunctions.writeOutfile(self.outfile, self.in_array, self.nams, self.removed)
        self.assertTrue(os.path.isfile(self.outfile))


class UtilityFunctionsCheckSeqTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    @parameterized.expand([
            ['./tests/test_files/example1.fasta', "nt"],
            ['./tests/test_files/example2.fasta', "aa"],
            ['./tests/test_files/uORF_nt.fasta', "nt"],
            ['./tests/test_files/uORF_aa.fasta', "aa"],
    ])
    def testSeqType(self, input, expected):
        in_array, names = readMSA(input)
        logger = logging.getLogger('path.to.module.under.test')
        type = utilityFunctions.seqType(in_array, logger)

        self.assertEqual(type, expected)

    @parameterized.expand([
            ['./tests/test_files/nonIUPAC.fasta']])
    def testSeqTypeNonIUPAC(self, input):
        arr, nams = utilityFunctions.FastaToArray(input)
        logger = logging.getLogger('path.to.module.under.test')
        with self.assertRaises(SystemExit) as cm:
            ali, nams = utilityFunctions.seqType(arr, logger)
        captured = self.capsys.readouterr()
        self.assertEqual(cm.exception.code, 1)
        self.assertTrue('Unknown nucleotides or amino acids detected'
                        in captured.out)



class UtilityFunctionsListFontsTest(unittest.TestCase):

    def setUp(self):
        self.outfile = "listFonts_test.png"

    def tearDown(self):
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

    def testListFonts(self):
        with warnings.catch_warnings():
            # Don't raise warnings for missing glyphs
            warnings.filterwarnings('ignore')
            plat = sys.platform
            if plat != "win32" and plat != "cygwin":
                utilityFunctions.listFonts(self.outfile)
                self.assertTrue(os.path.isfile(self.outfile))

class UtilityFunctionsUpdateNamsTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @parameterized.expand([
            [['Hi', 'I', 'Am', 'A', 'Test'], set(['Am']),
             ['Hi', 'I', 'A', 'Test']],
            [['Hi', 'I', 'Am', 'A', 'Test'], set(['The']),
             ['Hi', 'I', 'Am', 'A', 'Test']],
            [['Hi', 'I', 'Am', 'A', 'Test'],
             set(['Test', 'A', 'Am', 'I', 'Hi']),
             []]
            ])
    def testUpdateNams(self, nams, removed, expected):
        rm = utilityFunctions.updateNams(nams, removed)
        self.assertEqual(rm, expected)


class UtilityFunctionscheckArrLengthTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @parameterized.expand([[np.zeros([0, 0]), "removed all of the sequences"],
                           [np.zeros([0, 4]), "removed all of the sequences"],
                           [np.zeros([3]), "not all the same length"]])
    def testCheckArrLength0(self, arr, expected):
        logger = logging.getLogger('path.to.module.under.test')

        with self.assertRaises(SystemExit) as cm:
            utilityFunctions.checkArrLength(arr, logger)
        captured = self.capsys.readouterr()
        self.assertEqual(cm.exception.code, 1)
        self.assertTrue(expected in captured.out)

    @parameterized.expand([[np.zeros([10, 4])]])
    def testCheckArrLengthNot0(self, arr):
        logger = logging.getLogger('path.to.module.under.test')

        utilityFunctions.checkArrLength(arr, logger)
        captured = self.capsys.readouterr()

        self.assertTrue(len(captured.out) == 0)


class UtilityFunctionsRetainTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def setUp(self):
        self.arr, self.nams = utilityFunctions.FastaToArray(
            "tests/test_files/example1.fasta")

    def tearDown(self):
        pass

    @parameterized.expand([[None, [""], "", np.array([])],
                          [['Seq1'], [""], "", np.array(['Seq1'])],
                          [[''], ["S"], "", np.array(['Seq1',
                                                      'Seq2',
                                                      'Seq3',
                                                      'Seq4',
                                                      'Seq5',
                                                      'Seq6',
                                                      'Seq7'])],
                          [[''], [""], "tests/test_files/keep.txt",
                           np.array(['Seq2',
                                     'Seq3'])],
                          [['Seq1'], ["7"], "tests/test_files/keep.txt",
                           np.array(['Seq1', 'Seq2', 'Seq3', 'Seq7'])]])
    def testConfigRetainSeqs(self, retain, retainS, retainL, expected):
        logger = logging.getLogger('path.to.module.under.test')
        keeps = utilityFunctions.configRetainSeqs(retain, retainS, retainL,
                                                  self.nams, "func", logger,
                                                  False)
        self.assertTrue(np.array_equal(keeps, expected))

    def testConfigRetainSeqsFNF(self):
        logger = logging.getLogger('path.to.module.under.test')
        
        with self.assertRaises(FileNotFoundError):
            utilityFunctions.configRetainSeqs(None, [""], "hello.txt",
                                              self.nams, "func", logger,
                                              False)

    def testConfigRetainSeqsNoMatch(self):
        logger = logging.getLogger('path.to.module.under.test')

        utilityFunctions.configRetainSeqs(None, ["hello"], "",
                                          self.nams, "func", logger,
                                          False)
        captured = self.capsys.readouterr()
        self.assertTrue("No sequence names matching" in captured.out)

    @parameterized.expand([[['S1'], [""], ""],
                           [None, ["S"], "tests/test_files/keep_m.txt"]])
    def testConfigRetainSeqsNoMatchList(self, retain, retainS, retainL):
        logger = logging.getLogger('path.to.module.under.test')

        with self.assertRaises(RuntimeError):
            utilityFunctions.configRetainSeqs(retain, retainS, retainL,
                                              self.nams, "func", logger,
                                              False)


class UtilityFunctionsUpdateStartEndTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @parameterized.expand([[0, 10, [5, 6, 7], [0, 7]],
                          [10, 20, [], [10, 20]],
                          [10, 20, [5, 6, 15], [8, 17]]])
    def testConfigRetainSeqs(self, start, end, removed, expected):
        newstart, newend = utilityFunctions.updateStartEnd(start, end, removed)
        self.assertEqual(newstart, expected[0])
        self.assertEqual(newend, expected[1])
        