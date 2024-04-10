import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_utils.py")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_combine_class_ii_alleles(self):
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DQA1*06:02", "DQB1*06:43", "DQB1*06:44"])),
            sorted(["DQA1*06:02", "DQB1*06:43", "DQB1*06:44", "DQA1*06:02-DQB1*06:43", "DQA1*06:02-DQB1*06:44"])
        )

        #two different types of alleles
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DQA1*06:02", "DQB1*06:43", "DPA1*01:03", "DPB1*02:01"])),
            sorted(["DQA1*06:02", "DQB1*06:43", "DPA1*01:03", "DPB1*02:01", "DQA1*06:02-DQB1*06:43", "DPA1*01:03-DPB1*02:01"])
        )

        #DRA*01:01-DRB9*01:02 is invalid
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DRB9*01:02", "DRA*01:01"])),
            sorted(["DRB9*01:02", "DRA*01:01"])
        )
