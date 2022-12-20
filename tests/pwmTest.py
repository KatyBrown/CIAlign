#!/usr/bin/env python
'''
python3 -m unittest tests.pwmTest
in CIAlign folder
'''

import unittest
import logging
from parameterized import parameterized
import pandas as pd
import numpy as np
import os
import CIAlign.matrices as matrices
import CIAlign.utilityFunctions as utilityFunctions
import pandas.testing

class PWMTests(unittest.TestCase):

    def setUp(self):
        arr_rna, nams_rna = utilityFunctions.FastaToArray(
            "./tests/test_files/example1.fasta")
        arr_aa, nams_aa = utilityFunctions.FastaToArray(
            "./tests/test_files/example2.fasta")
        arr_dna, nams_dna = utilityFunctions.FastaToArray(
            "./tests/test_files/example3_section.fasta")

        self.arrs = {'rna': arr_rna,
                     'dna': arr_dna,
                     'aa': arr_aa}

        self.matrices = {'pfm': dict(), 'ppm': dict(), 'pwm': dict(),
                         'freq': dict(), 'alpha': dict()}

        for string in ['dna', 'rna', 'aa']:
            for mat in ['pfm', 'ppm', 'pwm']:
                tab = pd.read_csv("./tests/test_files/%s_%s.tsv" % (mat,
                                                                    string),
                                  sep="\t", header=None, index_col=0)
                tab.columns = pd.RangeIndex(np.shape(tab)[1])
                self.matrices[mat][string] = tab

            for stat in ['alpha', 'freq']:
                self.matrices[stat][string] = np.load(
                    "./tests/test_files/%s_%s.npy" % (stat, string))
    
        tab = pd.read_csv("./tests/test_files/pfm_dna_large.tsv",
                          sep="\t", header=None, index_col=0)
        tab.columns = pd.RangeIndex(np.shape(tab)[1])
        self.matrices['pfm']['dna_large'] = tab


    def tearDown(self):
        f = "./tests/test_files/test_meme.out"
        if os.path.exists(f):
            os.remove(f)


    @parameterized.expand([
        ['nt', ['A', 'C', 'G', 'T']],
        ['aa', ['A', 'C', 'D', 'E', 'F',
                'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R',
                'S', 'T', 'V', 'W', 'Y']]
         ])
    def testgetCoreRes(self, typ, expected):
        calculated = matrices.getCoreRes(typ)
        self.assertEqual(calculated, expected)


    @parameterized.expand([['nt', 'dna'],
                           ['aa', 'aa'],
                           ['nt', 'rna']])
    def testmakePFM(self, typ, subtyp):
        arr = self.arrs[subtyp]
        PFM, RNA = matrices.makePFM(arr, typ)
        pandas.testing.assert_frame_equal(PFM, self.matrices['pfm'][subtyp],
                                          check_names=False)


    @parameterized.expand([['nt', 'dna', False],
                           ['aa', 'aa', None],
                           ['nt', 'rna', True]])
    def testmakePPM(self, typ, subtyp, RNA):
        PFM = self.matrices['pfm'][subtyp]
        alpha = self.matrices['alpha'][subtyp]
        PPM = matrices.makePPM(PFM, alpha).round(3)
        pandas.testing.assert_frame_equal(PPM, self.matrices['ppm'][subtyp])


    @parameterized.expand([['nt', 'dna', False],
                           ['aa', 'aa', None],
                           ['nt', 'rna', True]])
    def testmakePWM(self, typ, subtyp, RNA):
        PPM = self.matrices['ppm'][subtyp]
        freq = self.matrices['freq'][subtyp]
        PWM = matrices.makePWM(PPM, freq)
        pandas.testing.assert_frame_equal(PWM, self.matrices['pwm'][subtyp])


    @parameterized.expand([['equal', 'nt', 'dna', False,
                            np.full((4, 48), 0.25)],
                           ['equal', 'aa', 'aa', None,
                            np.full((20, 1817), 0.05)],
                           ['equal', 'nt', 'rna', True, 
                            np.full((4, 96), 0.25)],
                           ['calc', 'nt', 'dna', False,
                           np.vstack([np.repeat(0.237, 48),
                                      np.repeat(0.142, 48),
                                      np.repeat(0.282, 48),
                                      np.repeat(0.338, 48)])],
                           ['calc2', 'nt', 'dna', True,
                            np.vstack([np.repeat(0.295, 48),
                                       np.repeat(0.169, 48),
                                       np.repeat(0.173, 48),
                                       np.repeat(0.364, 48)])],
                           ['xxx', 'nt', 'dna', True, np.zeros((1, 1))]])
    def testgetFreq(self, freqtype, typ, subtyp, RNA, expected):
        logger = logging.getLogger('path.to.module.under.test')
        PFM = self.matrices['pfm'][subtyp]
        if freqtype == 'calc2':
            PFM2 = self.matrices['pfm']['dna_large']
        else:
            PFM2 = None
        if freqtype == 'xxx':
            self.assertRaises(RuntimeError, matrices.getFreq, 
                              freqtype, logger, typ, RNA, PFM, PFM2)
        else:
            freqs = matrices.getFreq(freqtype, logger, typ, RNA, PFM,
                                     PFM2)
            self.assertTrue(np.array_equal(freqs, expected))


    @parameterized.expand([['calc', 'nt', 'dna', False, None,
                            np.load("./tests/test_files/alpha_dna_calc.npy")],
                           ['calc', 'aa', 'aa', None, None,
                            np.load("./tests/test_files/alpha_aa_calc.npy")],
                           ['calc', 'nt', 'rna', True, None,
                            np.load("./tests/test_files/alpha_rna_calc.npy")],
                           ['user', 'nt', 'dna', False, 1.0,
                            np.full((4, 48), 1.0)],
                           ['user', 'aa', 'aa', False, 1.5,
                            np.full((20, 1817), 1.5)],
                           ['xxx', 'aa', 'aa', False, 1.5,
                            np.zeros((1, 1))]])
    def testgetAlpha(self, alphatype, typ, subtyp, RNA, alphaval, expected):
        logger = logging.getLogger('path.to.module.under.test')
        PFM = self.matrices['pfm'][subtyp]
        freq = self.matrices['freq'][subtyp]
        if alphatype == 'xxx':
            self.assertRaises(RuntimeError, matrices.getAlpha, 
                              alphatype, logger, PFM, freq, alphaval)
        else:
            alpha = matrices.getAlpha(alphatype, logger, PFM, freq, alphaval)
            self.assertTrue(np.array_equal(alpha, expected))


    @parameterized.expand([['nt', 'dna', False,
                            "./tests/test_files/meme_dna.out"],
                           ['aa', 'aa', None,
                            "./tests/test_files/meme_aa.out"],
                           ['nt', 'rna', True,
                            "./tests/test_files/meme_rna.out"]])
    def testmemeFormat(self, typ, subtyp, RNA, expected):
        freq = self.matrices['freq'][subtyp]
        ppm = self.matrices['ppm'][subtyp]
        outfile = "./tests/test_files/test_meme.out"
        matrices.memeFormat(ppm, typ, RNA, freq, outfile, "")
        out = pd.read_csv(outfile, skiprows=8, sep=" ", header=None)
        exp = pd.read_csv(expected, skiprows=8, sep=" ", header=None)
        pandas.testing.assert_frame_equal(out, exp)