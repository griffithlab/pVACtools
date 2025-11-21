import unittest
import os
import re
import sys
import py_compile
import tempfile
from subprocess import PIPE
from subprocess import run as subprocess_run

from pvactools.tools.pvacseq import *
from tests.utils import *


def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'ml_predictor'
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
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq add_ml_predictions -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacseq',
            'add_ml_predictions.py'
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        input_file_I_aggregated = os.path.join(self.test_data_directory, 'MHC_Class_I', 'HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv')
        input_file_I_all_epitopes = os.path.join(self.test_data_directory, 'MHC_Class_I', 'HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv')
        input_file_II_aggregated = os.path.join(self.test_data_directory, 'MHC_Class_II', 'HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv')
        output_dir = tempfile.TemporaryDirectory()
        self.assertFalse(add_ml_predictions.main([
            '--class1-aggregated', input_file_I_aggregated,
            '--class1-all-epitopes', input_file_I_all_epitopes,
            '--class2-aggregated', input_file_II_aggregated,
            output_dir.name,
            'HCC1395'
        ]))