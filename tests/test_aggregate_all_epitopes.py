import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from lib.aggregate_all_epitopes import *
from .test_utils import *

class AggregateAllEptiopesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.executable    = os.path.join(cls.base_dir, "lib", "aggregate_all_epitopes.py")
        cls.test_data_dir = os.path.join(cls.base_dir, "tests", "test_data", "aggregate_all_epitopes")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_aggregate_all_epitopes_pvacseq_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(AggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.tsv"),
        ))

    def test_aggregate_all_epitopes_pvacfuse_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(AggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacfuse.tsv'), output_file.name, 'pVACfuse').execute())
        self.assertTrue(compare(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacfuse.tsv"),
        ))

    def test_aggregate_all_epitopes_pvacbind_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(AggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacbind.tsv'), output_file.name, 'pVACbind').execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacbind.tsv"),
        ))
