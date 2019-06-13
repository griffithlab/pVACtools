import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from lib.calculate_manufacturability import *

class CalculateManufacturabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.executable    = os.path.join(cls.base_dir, "lib", "calculate_manufacturability.py")
        cls.test_data_dir = os.path.join(cls.base_dir, "tests", "test_data", "calculate_manufacturability")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_calculate_manufacturability(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(CalculateManufacturability(os.path.join(self.test_data_dir, 'input.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.tsv"),
        ))
