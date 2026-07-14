import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools import *
from tests.utils import *

class PvacseqValidNetmhciipanVersionsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module_path = os.path.join(pvactools_directory(), "pvactools", "tools", "valid_netmhciipan_versions.py")

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
            'valid_netmhciipan_versions',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvactools valid_netmhciipan_versions -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_runs(self):
        self.assertFalse(valid_netmhciipan_versions.main([]))
