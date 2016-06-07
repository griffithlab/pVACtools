import unittest
import os
import sys
import py_compile
import subprocess

class pvacCompilationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

    def test_pvac_seq_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            "pVAC-Seq.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_main_compiles(self):
        compiled_main_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            "pvac_seq",
            "main.py"
        ))
        self.assertTrue(compiled_main_path)
