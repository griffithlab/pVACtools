import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from tempfile import TemporaryDirectory

from pvactools.tools.pvacsplice import *
from tests.utils import *

class PvacspliceDownloadExampleDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()

    def test_command(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacsplice',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'download_example_data',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacsplice download_example_data -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "download_example_data.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        self.assertFalse(download_example_data.main([output_dir.name]))
