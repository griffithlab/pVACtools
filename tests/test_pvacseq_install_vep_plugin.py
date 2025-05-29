import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from tempfile import TemporaryDirectory

from pvactools.tools.pvacseq import *
from tests.utils import *

class PvacseqInstallVepPluginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()

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
            'install_vep_plugin',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq install_vep_plugin -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "install_vep_plugin.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        self.assertFalse(install_vep_plugin.main([output_dir.name]))
