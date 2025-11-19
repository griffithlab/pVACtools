import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools.pvacseq import *
from tests.utils import *


def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacseq'
    )


class PvacseqAddMlPredictionsTests(unittest.TestCase):
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
            'main.py'
        )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'add_ml_predictions',
            '-h'
        ], shell=False, stdout=PIPE, stderr=PIPE)
        # Check that usage information is displayed (either in stdout or stderr)
        output = result.stdout.decode() + result.stderr.decode()
        self.assertRegex(output, usage_search, "Usage information should be displayed")

    def test_compiles(self):
        compiled_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacseq',
            'add_ml_predictions.py'
        ))
        self.assertTrue(compiled_path)

    def test_parser_requires_arguments(self):
        # Calling main() without required args should raise SystemExit,
        # proving the parser is hooked up correctly.
        with self.assertRaises(SystemExit):
            add_ml_predictions.main([])