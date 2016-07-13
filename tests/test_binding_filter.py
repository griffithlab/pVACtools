import unittest
import os
import tempfile
from filecmp import cmp
from subprocess import call
import sys
import py_compile

#python -m unittest tests/test_binding_filter.py
class BindingFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.binding_filter_path = os.path.join(cls.pVac_directory, "pvacseq", "lib", "binding_filter.py")
        cls.test_data_path= os.path.join(cls.pVac_directory, "tests", "test_data", "binding_filter")

    def test_binding_filter_runs_and_produces_expected_output(self):
        compiled_script_path = py_compile.compile(self.binding_filter_path)
        self.assertTrue(compiled_script_path)
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(call([
            sys.executable,
            compiled_script_path,
            os.path.join(
                self.test_data_path,
                'Test.HLA-A*29:02.9.netmhc.parsed.tsv'
            ),
            output_file.name
        ], shell=False))
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test_filtered.tsv"),
            False
        ))
