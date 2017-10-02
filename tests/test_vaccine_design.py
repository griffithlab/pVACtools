import unittest
import tempfile
import py_compile
import subprocess
import shutil
from filecmp import cmp
import os
import sys
import re
from subprocess import PIPE
from subprocess import run as subprocess_run

#python -m unittest tests/test_vaccine_design.py
class TestVaccineDesign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.python = sys.executable
        cls.executable = os.path.join(cls.base_dir, 'tools', 'pvacvector', 'run.py')
        cls.test_run_name = 'test_vaccine_design_produces_expected_output'
        cls.test_data_dir = os.path.join(cls.base_dir, 'tests', 'test_data', 'vaccine_design')
        cls.test_data_temp_dir = os.path.join(cls.test_data_dir, 'tmp')
        cls.input_tsv = os.path.join(cls.test_data_dir, 'input_parse_test_input.tsv')
        cls.input_vcf = os.path.join(cls.test_data_dir, 'input_parse_test_input.vcf')
        cls.input_file = os.path.join(cls.test_data_dir, 'Test.vaccine.results.input.fa')
        cls.method = 'ann'
        cls.keep_tmp = 'True'
        cls.allele = 'H-2-Kb'
        cls.epitope_length = '8'
        cls.seed = 'True'
        cls.input_n_mer = '25'

    def test_vaccine_design_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_pvacvector_compiles(self):
        compiled_path = py_compile.compile(os.path.join(
            self.base_dir,
            'tools',
            'pvacvector',
            'main.py'
        ))
        self.assertTrue(compiled_path)

    def test_pvacvectory_commands(self):
        pvac_script_path = os.path.join(
            self.base_dir,
            'tools',
            'pvacvector',
            'main.py'
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_vaccine_design_fa_input_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()

        call = subprocess.call([self.python,
                           self.executable,
                           self.test_run_name,
                           self.method,
                           self.allele,
                           '-f', self.input_file,
                           '-o', output_dir.name,
                           '-l', self.epitope_length,
                           '-n', self.input_n_mer,
                           '-k',
                           '-s'], shell=False)

        #vaccine design algorithm producing correct output with fasta input
        self.assertTrue(cmp(
            os.path.join(output_dir.name, self.test_run_name, self.test_run_name + '_results.fa'),
            os.path.join(self.test_data_dir, "Test.vaccine.results.output.fa")
        ))

        image_out = os.path.join(output_dir.name, self.test_run_name, 'vaccine.jpg')
        #vaccine visualization producing image
        self.assertTrue(os.path.exists(image_out))
        self.assertTrue(os.stat(image_out).st_size > 0)

        output_dir.cleanup()

    def test_vaccine_design_generate_fa_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()

        call = subprocess.call([self.python,
                           self.executable,
                           self.test_run_name,
                           self.method,
                           self.allele,
                           '-g',
                           '-t', self.input_tsv,
                           '-v', self.input_vcf,
                           '-o', output_dir.name,
                           '-l', self.epitope_length,
                           '-n', self.input_n_mer,
                           '-k',
                           '-s'], shell=False)

        #conversion from vcf to fasta file producing correct output, input file for vaccine design algorithm
        self.assertTrue(cmp(
                os.path.join(output_dir.name, self.test_run_name, "vaccine_design_input.fa"),
                os.path.join(self.test_data_dir, "input_parse_test_output.fa")
                ))

        image_out = os.path.join(output_dir.name, self.test_run_name, 'vaccine.jpg')
        #vaccine visualization producing image
        self.assertTrue(os.path.exists(image_out))
        self.assertTrue(os.stat(image_out).st_size > 0)

        output_dir.cleanup()
        
