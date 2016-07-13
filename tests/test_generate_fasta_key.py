import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class GenerateFastaKeyTests(unittest.TestCase):
    def setUp(self):
        self.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        self.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        self.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'generate_fasta_key')
        self.sample_name             = 'Test'
        self.peptide_sequence_length = 21

    def tearDown(self):
        del self.executable_dir
        del self.test_data_dir
        del self.sample_name
        del self.peptide_sequence_length

    def test_generate_fasta_key_runs_and_produces_expected_output(self):
        generate_fasta_key_input_file  = os.path.join(self.test_data_dir, ("%s_%s.fa" % (self.sample_name, self.peptide_sequence_length)))
        generate_fasta_key_output_file = tempfile.NamedTemporaryFile()
        generate_fasta_key_executable  = os.path.join(self.executable_dir, 'generate_fasta_key.py')
        self.assertTrue(py_compile.compile(generate_fasta_key_executable))

        self.assertFalse(call([
            self.python, generate_fasta_key_executable, generate_fasta_key_input_file, generate_fasta_key_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, ("%s_%s.key" % (self.sample_name, self.peptide_sequence_length)))
        self.assertTrue(cmp(generate_fasta_key_output_file.name, expected_output_file, shallow=False))

if __name__ == '__main__':
    unittest.main()
