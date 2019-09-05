#! /usr/bin/env python
import sys
sys.path.append("../src")
import glob
import CIAlign
import logging
import itertools
from collections import OrderedDict
import consensusSeq

class TestFunctions(object):
    def __init__(self):
        testfiles = glob.glob("test_alignments/test_cases/*fasta")
        realfiles = glob.glob("test_alignments/real_cases/*fasta")
        self.testfiles = testfiles
        self.realfiles = realfiles
        self.files = testfiles + realfiles
        log = logging.getLogger(__name__)
        log.setLevel(logging.INFO)
    
        handler = logging.FileHandler("test_alignments/log.log")
        handler.setLevel(logging.INFO)
    
        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.log = log

    def run_test(self, params_test, function_test, testname, direc='cols'):
        for i, file in enumerate(self.files):
            print (i, testname, file)
            pvals = list(params_test.values())
            for p in itertools.product(*pvals):
                print (*p)
                arr, nams = CIAlign.FastaToArray(file)
                if direc == 'cols':
                    positions = list(range(0, len(arr[0])))
                    arr2, r2 = function_test(arr, positions, rmfile, self.log, *p) 
                elif direc == 'rows':
                    arr2, r2 = function_test(arr, nams, rmfile, self.log, *p) 

            
            
    def test_remove_insertions(self):
        params_test = OrderedDict()
        params_test['insertion_min_size'] = [1, 10, 50]
        params_test['insertion_max_size'] = [1, 10, 50]
        params_test['insertion_min_flank'] = [5, 25]

        self.run_test(params_test, parsingFunctions.removeInsertions, 'remove_insertions',
                      direc='cols')        
        
    def test_remove_short(self):
        params_test = OrderedDict()
        params_test['remove_min_length'] = [10, 50, 100]
        self.run_test(params_test, parsingFunctions.removeTooShort, 'remove_short', direc='rows')


    def test_make_consensus(self):
        params_test = OrderedDict()
        params_test['consensus_type'] = ['majority', 'majority_nongap']
        self.run_test(params_test, parsingFunctions.findConsensus, 'make_consensus',
                      direc='NA')

        
    def test_crop_ends(self):
        params_test = OrderedDict()
        params_test['crop_ends_mingap'] = [0, 5, 10, 20, 50]
        self.run_test(params_test, parsingFunctions.cropEnds, 'crop_ends',
                      direc='rows')        
  

    def test_remove_badlyaligned(self):
        params_test = OrderedDict()
        params_test['remove_badlyaligned_minperc'] = [0, 0.5, 0.8, 1.0]
        self.run_test(params_test, parsingFunctions.removeBadlyAligned, 'remove_badlyaligned', direc='rows')
    
    """
    def test_remove_gaponly(self):
        pass
        
    def test_make_similarity_matrix_input(self):
        params_test = ['make_simmatrix_dp', 'make_simmatrix_minoverlap',
                       'make_simmatrix_keepgaps']

    def test_make_similarity_matrix_output(self):
        params_test = ['make_simmatrix_dp', 'make_simmatrix_minoverlap',
                       'make_simmatrix_keepgaps']
    
    def test_plot_input(self):
        params_test = ['plot_dpi', 'plot_format', 'plot_width', 'plot_height']

    def test_plot_output(self):
        params_test = ['plot_dpi', 'plot_format', 'plot_width', 'plot_height']

    def test_plot_markup(self):
        params_test = ['plot_dpi', 'plot_format', 'plot_width', 'plot_height']
        
    def test_make_sequence_logo(self):
        params_test = ['sequence_logo_type', 'sequence_logo_dpi',
                       'sequence_logo_font']
"""        
def main():
    t = TestFunctions()
    # t.test_remove_insertions()
    # t.test_remove_short()
    # t.test_make_consensus()
    t.test_crop_ends()


if __name__ == "__main__":
    main()
