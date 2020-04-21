import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from lib.rank_epitopes import *

class RankEpitopesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.executable    = os.path.join(cls.base_dir, "lib", "rank_epitopes.py")
        cls.test_data_dir = os.path.join(cls.base_dir, "tests", "test_data", "rank_epitopes")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_rank_epitopes_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(RankEpitopes(os.path.join(self.test_data_dir, 'input.tsv'), output_file.name, 'median').execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.tsv"),
        ))
