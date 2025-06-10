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

class PvacseqNetChopTests(unittest.TestCase):
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
            'net_chop',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq net_chop -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "net_chop.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        input_file = os.path.join(self.test_data_directory, 'Test.all_epitopes.short.tsv')
        input_fasta = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.fasta')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(net_chop.main([
            input_file,
            input_fasta,
            output_file.name,
        ]))
