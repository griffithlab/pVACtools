import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools.pvacfuse import *
from tests.utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacfuse'
    )

class PvacfuseUpdateTiersTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()

    def test_command(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacfuse',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'update_tiers',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacfuse update_tiers -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacfuse",
            "update_tiers.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        input_file = os.path.join(self.test_data_directory, 'arriba_fusions', 'MHC_Class_I', 'Test.all_epitopes.aggregated.tsv')
        self.assertFalse(update_tiers.main([input_file]))
