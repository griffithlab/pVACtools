import unittest
import os
import tempfile
from filecmp import cmp
from subprocess import call
import sys
import py_compile

#python -m unittest tests/test_coverage_filter.py
class CoverageFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.executable    = os.path.join(cls.base_dir, "pvacseq", "lib", "coverage_filter.py")
        cls.test_data_dir = os.path.join(cls.base_dir, "tests", "test_data", "coverage_filter")

    def test_binding_filter_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        coverage_filter_cmd = "%s  %s  %s %s" % (
            self.python,
            self.executable,
            os.path.join(
                self.test_data_dir,
                'Test_filtered.readcounts.tsv'
            ),
            output_file.name
        )
        self.assertFalse(call([coverage_filter_cmd], shell=True))
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test_filtered.readcounts.covfilt.tsv"),
        ))
