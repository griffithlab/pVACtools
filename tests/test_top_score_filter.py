import unittest
import sys
import os
import py_compile
import tempfile
from filecmp import cmp

from pvactools.lib.top_score_filter import PvacseqTopScoreFilter, PvacfuseTopScoreFilter, PvacbindTopScoreFilter, PvacspliceTopScoreFilter
from tests.utils import *

class TopScoreFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'top_score_filter.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'top_score_filter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_runs_and_creates_expected_file_median(self):
        input_file = os.path.join(self.test_data_dir, 'input.tsv')
        output_file = tempfile.NamedTemporaryFile()

        PvacseqTopScoreFilter(input_file, output_file.name).execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_median.tsv')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_lowest(self):
        input_file = os.path.join(self.test_data_dir, 'input.tsv')
        output_file = tempfile.NamedTemporaryFile()

        PvacseqTopScoreFilter(input_file, output_file.name, top_score_metric='lowest').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_lowest.tsv')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_fusion(self):
        input_file = os.path.join(self.test_data_dir, 'input_fusion.tsv')
        output_file = tempfile.NamedTemporaryFile()

        PvacfuseTopScoreFilter(input_file, output_file.name, top_score_metric='median').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_fusion.tsv')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_pvacbind(self):
        input_file = os.path.join(self.test_data_dir, 'input_pvacbind.tsv')
        output_file = tempfile.NamedTemporaryFile()

        PvacbindTopScoreFilter(input_file, output_file.name, top_score_metric='median').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_pvacbind.tsv')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_pvacsplice(self):
        input_file = os.path.join(self.test_data_dir, 'input_pvacsplice.tsv')
        output_file = tempfile.NamedTemporaryFile()

        PvacspliceTopScoreFilter(input_file, output_file.name, top_score_metric='median').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_pvacsplice.tsv')
        self.assertTrue(cmp(output_file.name, expected_output_file))
