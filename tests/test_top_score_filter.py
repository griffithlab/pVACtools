import unittest
import sys
import os
import py_compile
import tempfile
from lib.top_score_filter import *

def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result

class TopScoreFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'top_score_filter.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'top_score_filter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_runs_and_creates_expected_file_median(self):
        input_file = os.path.join(self.test_data_dir, 'input.tsv')
        output_file = tempfile.NamedTemporaryFile()

        TopScoreFilter(input_file, output_file.name, 'median').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_median.tsv')
        self.assertTrue(compare(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_lowest(self):
        input_file = os.path.join(self.test_data_dir, 'input.tsv')
        output_file = tempfile.NamedTemporaryFile()

        TopScoreFilter(input_file, output_file.name, 'lowest').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_lowest.tsv')
        self.assertTrue(compare(output_file.name, expected_output_file))

    def test_runs_and_creates_expected_file_fusion(self):
        input_file = os.path.join(self.test_data_dir, 'input_fusion.tsv')
        output_file = tempfile.NamedTemporaryFile()

        TopScoreFilter(input_file, output_file.name, 'median').execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_fusion.tsv')
        self.assertTrue(compare(output_file.name, expected_output_file))
