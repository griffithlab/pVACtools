import unittest
import os
import sys
import tempfile
from filecmp import cmp
import py_compile
import pandas as pd

from pvactools.lib.filter_regtools_results import FilterRegtoolsResults
from tests.utils import *

#python -m unittest tests/test_pvacsplice_filter_regtools_results.py
class FilterRegtoolsResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "filter_regtools_results.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")
        # inputs

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_filter_regtools_results_runs_and_produces_expected_output(self):
        input_file = os.path.join(self.test_data_dir, 'inputs', 'splice_junctions_chr1.tsv')
        output_dir = tempfile.TemporaryDirectory()
        gtf_df = pd.read_csv(os.path.join(self.test_data_dir, 'results', 'Test.gtf.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str})

        # iterate free params
        for score in [5, 10, 50]:
            for distance in [25, 50, 100]:
                output_file = os.path.join(output_dir.name, 'sample_{}_{}_filtered.tsv'.format(score, distance))
                params = {
                    'input_file'  : input_file,
                    'output_file' : output_file,
                    'gtf_data'    : gtf_df,
                    'score'       : score,
                    'distance'    : distance,
                }
                FilterRegtoolsResults(**params).execute()

                expected_file = os.path.join(
                    self.test_data_dir,
                    'results',
                    'Test.{}_{}_filtered.tsv'.format(score, distance)
                )
                # direct file comparison to existing
                self.assertTrue(cmp(
                        output_file,
                        expected_file,
                        False
                    ),
                    "files don't match {} - {}".format(output_file, expected_file)
                )

        output_dir.cleanup()
