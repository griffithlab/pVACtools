import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools import *
from tests.utils import *

class PvacseqAlleleSpecificCutoffsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module_path = os.path.join(pvactools_directory(), "pvactools", "tools", "allele_specific_cutoffs.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.module_path))

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
            'allele_specific_cutoffs',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvactools allele_specific_cutoffs -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_runs(self):
        self.assertFalse(allele_specific_cutoffs.main([]))
