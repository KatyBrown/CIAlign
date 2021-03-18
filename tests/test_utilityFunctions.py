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

# class utilityFunctions


class utilityFunctionsgetColoursTests(unittest.TestCase):

    def test_getAAColours(self):

        AAcolours = utilityFunctions.getAAColours()

        self.assertEqual(len(AAcolours), 26)

    def test_getNtColours(self):

        Ntcolours = utilityFunctions.getNtColours()

        self.assertEqual(len(Ntcolours), 18)


class utilityFunctionsMSAInputTests(unittest.TestCase):

    def setUp(self):

        self.input = "./tests/test_files/example1.fasta"

        self.in_array = []
        self.nams = []
        input_ali = AlignIO.read(open(self.input), "fasta")
        for record in input_ali:
            self.in_array.append(record.seq)
            self.nams.append(record.id)
        self.in_array = np.array(self.in_array)


    def tearDown(self):
        pass

    def test_replaceUbyT(self):

        result_ali = utilityFunctions.replaceUbyT(self.in_array)
        findU = np.where(result_ali == "U")
        findu = np.where(result_ali == "u")

        self.assertFalse(len(findU[0]) > 0)
        self.assertFalse(len(findU[1]) > 0)
        self.assertFalse(len(findu[0]) > 0)
        self.assertFalse(len(findu[1]) > 0)

    def test_unAlign(self):

        result_ali = utilityFunctions.unAlign(self.in_array)
        findGap = np.where(result_ali == "-")

        self.assertFalse(len(findGap[0]) > 0)
        self.assertFalse(len(findGap[1]) > 0)

    def test_FastaToArray(self):

        logger = logging.getLogger('path.to.module.under.test')
        with mock.patch.object(logger, 'debug') as mock_debug:
            ali, nams = utilityFunctions.FastaToArray(self.input, logger)

        # self.assertEqual(nams.size, self.nams.size)
        self.assertEqual(ali[0,:].size, self.in_array[0,:].size)
        self.assertEqual(len(self.in_array), len(ali))
        self.assertEqual(len(nams), len(self.nams))
        self.assertTrue((ali == self.in_array).all())
        self.assertTrue(nams == self.nams)


class ultilityFunctionsWriteOutfileTest(unittest.TestCase):

    def setUp(self):

        self.input = "./tests/test_files/example1.fasta"

        self.in_array = []
        self.nams = []
        input_ali = AlignIO.read(open(self.input), "fasta")
        for record in input_ali:
            self.in_array.append(record.seq)
            self.nams.append(record.id)
        self.in_array = np.array(self.in_array)


    def tearDown(self):
        pass

    def test_writeOutfile(self):
