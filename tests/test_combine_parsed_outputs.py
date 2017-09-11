import unittest
import sys
import os
import tempfile
import py_compile
from subprocess import call
from filecmp import cmp

class CombineParsedOutputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python         = sys.executable
        base_dir           = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'combine_parsed_outputs.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'combine_parsed_outputs')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_combine_parsed_outputs_generates_expected_files(self):
        combine_parsed_outputs_output_file = tempfile.NamedTemporaryFile()
        combine_parsed_outputs_command = [
            self.python,
            self.executable,
            os.path.join(self.test_data_dir, 'Test.HLA-E*01:01.9.parsed.tsv'),
            os.path.join(self.test_data_dir, 'Test.HLA-G*01:09.9.parsed.tsv'),
            combine_parsed_outputs_output_file.name,
        ]
        self.assertFalse(call(combine_parsed_outputs_command))

        expected_output_file  = os.path.join(self.test_data_dir, "Test.combined.parsed.tsv")
        self.assertTrue(cmp(combine_parsed_outputs_output_file.name, expected_output_file))
