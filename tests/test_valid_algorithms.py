import unittest
import os
import sys
import py_compile
from unittest.mock import patch
import io

from pvactools.lib.valid_algorithms import ValidAlgorithms
from tests.utils import *

class ValidAlgorithmsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.valid_algorithms_path = os.path.join(pvactools_directory(), "pvactools", "lib", "valid_algorithms.py")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.valid_algorithms_path))

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_allele(self, mock_stdout):
        self.assertFalse(ValidAlgorithms(
            "HLA-A*01:01",
            None
        ).print_valid_algorithms())
        self.assertNotIn("NetMHCIIpan", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_species(self, mock_stdout):
        self.assertFalse(ValidAlgorithms(
            None,
            "mouse"
        ).print_valid_algorithms())
        self.assertNotIn("BigMHC_IM", mock_stdout.getvalue())

    def test_allele_species_mismatch(self):
        with self.assertRaises(Exception) as context:
            self.assertFalse(ValidAlgorithms(
                "HLA-A*02:01",
                "mouse"
            ).print_valid_algorithms())
        self.assertIn("Given species does not match given allele.", str(context.exception))

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_valid_algorithms_with_alleles_and_species(self, mock_stdout):
        self.assertFalse(ValidAlgorithms(
            "HLA-A*02:01",
            "human"
        ).print_valid_algorithms())
        self.assertIn("NetMHC", mock_stdout.getvalue())
