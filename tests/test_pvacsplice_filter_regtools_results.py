import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

from pvactools.tools.pvacsplice.filter_regtools_results import FilterRegtoolsResults
from .test_utils import *

class FilterRegtoolsResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "filter_regtools_results.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_filter_regtools_results_runs_and_produces_expected_output(self):
        input_file = os.path.join(self.test_data_dir, 'Test_junctions.tsv')
        output_dir = tempfile.TemporaryDirectory()

        # iterate free params
        for score in [10, 100, 1000]:
            for dist in [50, 100, 150]:
                output_file   = os.path.join(output_dir.name, 'sample_{}_{}_filtered.tsv'.format(score, dist))
                #output_file   = os.path.join(self.test_data_dir, "filter_regtools_results", 'Test_{}_{}_filtered.tsv'.format(score, dist))
                #do below and perhaps manually check scores & distances within expected range?
                params = {
                    'input_file'  : input_file,
                    'output_file' : output_file,
                    'score'       : score,
                    'distance'    : dist,
                }
                filtered = FilterRegtoolsResults(**params)
                filtered.execute()

                #direct file comparison to existing
                expected_file = os.path.join(self.test_data_dir, "filter_regtools_results", 'Test_{}_{}_filtered.tsv'.format(score, dist))
                self.assertTrue(compare(output_file, expected_file), "files don't match {} - {}".format(output_file, expected_file))

        output_dir.cleanup()
