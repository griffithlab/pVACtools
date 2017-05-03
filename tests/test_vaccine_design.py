import unittest
import tempfile
import py_compile
import subprocess
import shutil
from filecmp import cmp
import os
import sys

#python -m unittest tests/test_vaccine_design.py
class TestVaccineDesign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.python = sys.executable
        cls.executable_dir = os.path.join(cls.base_dir, 'pvacseq', 'lib')
        cls.executable = os.path.join(cls.executable_dir, 'vaccine_design.py')
        cls.test_run_name = 'test_vaccine_design_produces_expected_output'
        cls.test_data_dir = os.path.join(cls.base_dir, 'tests', 'test_data', 'vaccine_design')
        cls.test_data_temp_dir = os.path.join(cls.test_data_dir, 'tmp')
        cls.input_file = os.path.join(cls.test_data_dir, 'Test.vaccine.results.input.fa')
        cls.method = 'ann'
        cls.keep_tmp = 'True'
        cls.allele = 'H-2-Kb'
        cls.epitope_length = '8'
        cls.seed = 'True'

    def test_vaccine_design_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_vaccine_design_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()

        call = subprocess.call([self.python,
                           self.executable,
                           self.test_run_name,
                           self.input_file,
                           self.method,
                           self.allele,
                           '-o', output_dir.name,
                           '-l', self.epitope_length,
                           '-k',
                           '-s'], shell=False)

        self.assertTrue(cmp(
            os.path.join(output_dir.name, self.test_run_name, self.test_run_name + '_results.fa'),
            os.path.join(self.test_data_dir, "Test.vaccine.results.output.fa")
        ))

        output_dir.cleanup()
        
