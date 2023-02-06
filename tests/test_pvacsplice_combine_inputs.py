import unittest
import os
import sys
import tempfile
from filecmp import cmp
import pandas as pd
import py_compile

from pvactools.lib.combine_inputs import CombineInputs
from tests.utils import *

class CombineInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "combine_inputs.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice", 'results')
        # inputs 
        cls.junctions_df = pd.read_csv(os.path.join(cls.test_data_dir, 'Test.10_100_filtered.tsv'), sep='\t') # default filters
        cls.variant_file = os.path.join(cls.test_data_dir, 'Test.annotated.tsv')

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_combine_inputs_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_combined.tsv')
        params = {
            'junctions_df' : self.junctions_df,
            'variant_file' : self.variant_file,
            'output_file'  : output_file,
            'output_dir'   : output_dir
        }
        CombineInputs(**params).execute()
        expected_file = os.path.join(self.test_data_dir, 'Test.combined.tsv')
        self.assertTrue(cmp(
            output_file, 
            expected_file), 
            "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()
