import unittest
import os
import sys
import py_compile
from unittest.mock import patch
import io
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools import *
from tests.utils import *

class ValidAlgorithmsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.valid_algorithms_path = os.path.join(pvactools_directory(), "pvactools", "tools", "valid_algorithms.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.valid_algorithms_path))

    def test_command(self):
        pvac_script_path = os.path.join(
            pvactools_directory(),
            'pvactools',
            'tools',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'valid_algorithms',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvactools valid_algorithms -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_allele(self, mock_stdout):
        self.assertFalse(valid_algorithms.main([
            "--allele",
            "HLA-A*01:01",
        ]))
        self.assertNotIn("NetMHCIIpan", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_species(self, mock_stdout):
        self.assertFalse(valid_algorithms.main([
            "--species",
            "mouse"
        ]))
        self.assertNotIn("BigMHC_IM", mock_stdout.getvalue())

    def test_allele_species_mismatch(self):
        with self.assertRaises(Exception) as context:
            self.assertFalse(valid_algorithms.main([
                "--allele",
                "HLA-A*02:01",
                "--species",
                "mouse"
            ]))
        self.assertIn("Given species does not match given allele.", str(context.exception))

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_alleles_and_species(self, mock_stdout):
        self.assertFalse(valid_algorithms.main([
            "--allele",
            "HLA-A*02:01",
            "--species",
            "human"
        ]))
        self.assertIn("NetMHC", mock_stdout.getvalue())
