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

class ValidAllelesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.valid_alleles_path = os.path.join(pvactools_directory(), "pvactools", "tools", "valid_alleles.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.valid_alleles_path))

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
            'valid_alleles',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvactools valid_alleles -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_algorithm(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--prediction-algorithm",
            "NetMHCpan",
            "--species",
            "human",
        ]))
        self.assertNotIn("DPA", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_species(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--species",
            "mouse",
        ]))
        self.assertNotIn("HLA", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_algorithm_and_species(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--prediction-algorithm",
            "NetMHC",
            "--species",
            "mouse",
        ]))
        self.assertIn("H-2-Db", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_length(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--species",
            "human",
            "--length",
            "8"
        ]))
        self.assertIn("HLA", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_algorithm_and_length(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--prediction-algorithm",
            "NetMHCpan",
            "--species",
            "human",
            "--length",
            "8"
        ]))
        self.assertIn("HLA", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_species_and_length(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--species",
            "mouse",
            "--length",
            "8"
        ]))
        self.assertIn("H-2-Db", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_alleles_with_algorithm_and_species_and_length(self, mock_stdout):
        self.assertFalse(valid_alleles.main([
            "--prediction-algorithm",
            "NetMHC",
            "--species",
            "mouse",
            "--length",
            "8"
        ]))
        self.assertIn("H-2-Db", mock_stdout.getvalue())
