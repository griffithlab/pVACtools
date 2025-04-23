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

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_combine_inputs_runs_and_produces_expected_output(self):
        junctions_df = pd.read_csv(os.path.join(self.test_data_dir, 'Test.10_100_filtered.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str}) # default filters
        variant_file = os.path.join(self.test_data_dir, 'Test.annotated.tsv')
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_combined.tsv')
        params = {
            'junctions_df' : junctions_df,
            'variant_file' : variant_file,
            'output_file'  : output_file,
            'output_dir'   : output_dir.name,
        }
        combined_df = CombineInputs(**params).execute()
        combined_df.to_csv(output_file, sep='\t', index=False)

        expected_file = os.path.join(self.test_data_dir, 'Test.combined.tsv')
        self.assertTrue(cmp(
                output_file,
                expected_file),
                "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()

    def test_combine_inputs_GENCODE_runs_and_produces_expected_output(self):
        junctions_df = pd.read_csv(os.path.join(self.test_data_dir, 'Test.GENCODE.5_100_filtered.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str}) # default filters
        variant_file = os.path.join(self.test_data_dir, 'Test.GENCODE.annotated.tsv')
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_combined.GENCODE.tsv')
        params = {
            'junctions_df' : junctions_df,
            'variant_file' : variant_file,
            'output_file'  : output_file,
            'output_dir'   : output_dir.name,
        }
        combined_df = CombineInputs(**params).execute()
        combined_df.to_csv(output_file, sep='\t', index=False)

        expected_file = os.path.join(self.test_data_dir, 'Test.combined.GENCODE.tsv')
        self.assertTrue(cmp(
                output_file,
                expected_file),
                "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()

    def test_combine_inputs_runs_and_produces_expected_output_mouse(self):
        mouse_junctions_df = pd.read_csv(os.path.join(self.test_data_dir, 'results_mouse', 'Test.filtered.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str}) # default filters
        mouse_variant_file = os.path.join(self.test_data_dir, 'results_mouse', 'Test.annotated.tsv')
        output_dir = tempfile.TemporaryDirectory() 
        output_file = os.path.join(output_dir.name, 'sample_combined.tsv')
        params = {
            'junctions_df' : mouse_junctions_df,
            'variant_file' : mouse_variant_file,
            'output_file'  : output_file,
            'output_dir'   : output_dir.name,
        }
        combined_df = CombineInputs(**params).execute()
        combined_df.to_csv(output_file, sep='\t', index=False)

        expected_file = os.path.join(self.test_data_dir,'results_mouse', 'Test.combined.tsv')
        self.assertTrue(cmp(
                output_file,
                expected_file),
                "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()
