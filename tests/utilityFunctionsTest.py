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

import CIAlign.utilityFunctions as utilityFunctions
from tests.helperFunctions import readMSA


class UtilityFunctionsgetColoursTests(unittest.TestCase):

    def testGetAAColours(self):
        AAcolours = utilityFunctions.getAAColours()

        self.assertEqual(len(AAcolours), 26)

    def testGetNtColours(self):
        Ntcolours = utilityFunctions.getNtColours()

        self.assertEqual(len(Ntcolours), 18)


class UtilityFunctionsMSAInputTests(unittest.TestCase):

    def setUp(self):
        self.input = "./tests/test_files/example1.fasta"
        self.in_array, self.nams = readMSA(self.input)

    def tearDown(self):
        pass

    def testReplaceUbyT(self):
        result_ali = utilityFunctions.replaceUbyT(self.in_array)
        findU = np.where(result_ali == "U")
        findu = np.where(result_ali == "u")

        self.assertFalse(len(findU[0]) > 0)
        self.assertFalse(len(findU[1]) > 0)
        self.assertFalse(len(findu[0]) > 0)
        self.assertFalse(len(findu[1]) > 0)

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

        type = utilityFunctions.seqType(in_array)

        self.assertEqual(type, expected)

class UtilityFunctionsListFontsTest(unittest.TestCase):

    def setUp(self):
        self.outfile = "listFonts_test.png"

    def tearDown(self):
        os.remove(self.outfile)

    def testListFonts(self):
        utilityFunctions.listFonts(self.outfile)

        self.assertTrue(os.path.isfile(self.outfile))
