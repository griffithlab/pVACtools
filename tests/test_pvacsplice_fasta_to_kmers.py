import unittest
import os
import sys
import tempfile
from filecmp import cmp
import py_compile
import pandas as pd

from pvactools.lib.fasta_to_kmers import *
from tests.utils import *

#python -m unittest tests/test_pvacsplice_filter_regtools_results.py
class FastaToKmersTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "fasta_to_kmers.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_fasta_to_kmers_runs_and_produces_expected_output(self):
        tscript_fasta = os.path.join(self.test_data_dir, 'results', 'Test.transcripts.fa')
        class_i_epitope_length = [8,9,10,11]
        class_ii_epitope_length = [12,13,14,15,16]
        output_dir = tempfile.TemporaryDirectory()
        #output_file = os.path.join(output_dir.name, f'sample.{l}.kmers.tsv')
        params = {
            'fasta'           : tscript_fasta,
            'output_dir'      : output_dir.name,
            'class_i_epitope_length' : class_i_epitope_length,
            'class_ii_epitope_length': class_ii_epitope_length,
            'class_i_hla'     : 'HLA-A',
            'class_ii_hla'    : 'DBP1',
            'sample_name'     : 'sample',
        }           
        fasta = FastaToKmers(**params)
        fasta.execute()

        for x in class_i_epitope_length + class_ii_epitope_length:
            expected_file = os.path.join(self.test_data_dir, 'results', f'Test.{x}_kmers.fa')
            output_file = os.path.join(output_dir.name, f'sample.{x}.fa')

            self.assertTrue(cmp(
                    output_file, 
                    expected_file, 
                    False
                ), 
                "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()