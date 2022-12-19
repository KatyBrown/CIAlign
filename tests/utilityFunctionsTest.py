#! /usr/bin/env python
import logging
import unittest
from unittest import mock
from mock import patch
from parameterized import parameterized, parameterized_class
import matplotlib.font_manager
import sys
import logging
import numpy as np
from Bio import AlignIO
import os
from os import path
import sys
import io
import warnings

import CIAlign

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
        result_ali = utilityFunctions.replaceUbyT(self.in_array, rev=False)
        findU = np.where(result_ali == "U")
        findu = np.where(result_ali == "u")

        self.assertFalse(len(findU[0]) > 0)
        self.assertFalse(len(findU[1]) > 0)
        self.assertFalse(len(findu[0]) > 0)
        self.assertFalse(len(findu[1]) > 0)

    def testReplaceTbyU(self):
        result_ali = utilityFunctions.replaceUbyT(self.in_array, rev=True)
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

    def setUp(self):
        self.input = "./tests/test_files/empty.fasta"

    def tearDown(self):
        pass

    @parameterized.expand([
            ['./tests/test_files/corrupt.fasta']])
    def testFastaToArrayCorrupt(self, ali_path):
        logger = logging.getLogger('path.to.module.under.test')

        # capture output code from @paxdiablo on
        # https://stackoverflow.com/questions/33767627/python-write-unittest-for-console-print
        capturedOutput = io.StringIO()
        capturedError = io.StringIO()
        sys.stdout = capturedOutput
        sys.stderr = capturedError
        with self.assertRaises(SystemExit) as cm:
            ali, nams = utilityFunctions.FastaToArray(ali_path, logger)
        self.assertEqual(cm.exception.code, 1)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        printed_stdout = capturedOutput.getvalue()
        printed_stderr = capturedError.getvalue()
        
        self.assertTrue('needs to be in FASTA format' in printed_stdout)
        self.assertTrue('needs to be in FASTA format' in printed_stderr)


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

        capturedOutput = io.StringIO()
        capturedError = io.StringIO()
        sys.stdout = capturedOutput
        sys.stderr = capturedError
        with self.assertRaises(SystemExit) as cm:
            ali, nams = utilityFunctions.seqType(arr, logger)
        self.assertEqual(cm.exception.code, 1)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        printed_stdout = capturedOutput.getvalue()
        printed_stderr = capturedError.getvalue()
        
        self.assertTrue('Unknown nucleotides or amino acids detected'
                        in printed_stdout)
        self.assertTrue('Unknown nucleotides or amino acids detected'
                        in printed_stderr)


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

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @parameterized.expand([[np.zeros([0, 0]), "removed all of the sequences"],
                           [np.zeros([0, 4]), "removed all of the sequences"],
                           [np.zeros([3]), "not all the same length"]])
    def testCheckArrLength0(self, arr, expected):
        logger = logging.getLogger('path.to.module.under.test')
        # capture output code from @paxdiablo on
        # https://stackoverflow.com/questions/33767627/python-write-unittest-for-console-print
        capturedOutput = io.StringIO()
        capturedError = io.StringIO()
        sys.stdout = capturedOutput
        sys.stderr = capturedError
        with self.assertRaises(SystemExit) as cm:
            utilityFunctions.checkArrLength(arr, logger)
        self.assertEqual(cm.exception.code, 1)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        printed_stdout = capturedOutput.getvalue()
        printed_stderr = capturedError.getvalue()
        
        self.assertTrue(expected in printed_stdout)
        self.assertTrue(expected in printed_stderr)

    @parameterized.expand([[np.zeros([10, 4])]])
    def testCheckArrLengthNot0(self, arr):
        logger = logging.getLogger('path.to.module.under.test')
        # capture output code from @paxdiablo on
        # https://stackoverflow.com/questions/33767627/python-write-unittest-for-console-print
        capturedOutput = io.StringIO()
        capturedError = io.StringIO()
        sys.stdout = capturedOutput
        sys.stderr = capturedError

        utilityFunctions.checkArrLength(arr, logger)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        printed_stdout = capturedOutput.getvalue()
        printed_stderr = capturedError.getvalue()
        
        self.assertTrue(len(printed_stdout) == 0)
        self.assertTrue(len(printed_stderr) == 0)


class UtilityFunctionsRetainTest(unittest.TestCase):

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
    
        capturedOutput = io.StringIO()
        capturedError = io.StringIO()
        sys.stdout = capturedOutput
        sys.stderr = capturedError

        utilityFunctions.configRetainSeqs(None, ["hello"], "",
                                          self.nams, "func", logger,
                                          False)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        printed_stdout = capturedOutput.getvalue()
        printed_stderr = capturedError.getvalue()
        
        self.assertTrue("No sequence names matching" in printed_stdout)
        self.assertTrue("No sequence names matching" in printed_stderr)

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
        