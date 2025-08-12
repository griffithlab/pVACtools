import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from tempfile import NamedTemporaryFile

from pvactools.tools.pvacseq import *
from tests.utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacseq'
    )

class PvacseqIdentifyProblematicAminoAcidsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()

    def test_command(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacseq',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'identify_problematic_amino_acids',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq identify_problematic_amino_acids -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "identify_problematic_amino_acids.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.MHC_I.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(identify_problematic_amino_acids.main([input_file, output_file.name, "C"]))
