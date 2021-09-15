import unittest
import os
import tempfile
from filecmp import cmp
from subprocess import call
import sys
import py_compile

from .test_utils import *

#python -m unittest tests/test_coverage_filter.py
class CoverageFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacseq", "coverage_filter.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "coverage_filter")

    def test_coverage_filter_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(call([
            self.python,
            self.executable,
            os.path.join(
                self.test_data_dir,
                'Test.combined.parsed.tsv'
            ),
            output_file.name
        ], shell=False))
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test.filtered.coverage.tsv"),
        ))
