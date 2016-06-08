import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class GenerateFastaTests(unittest.TestCase):
    def setUp(self):
        self.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        self.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        self.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'generate_fasta')
        self.sample_name             = 'Test'
        self.peptide_sequence_length = 21

    def tearDown(self):
        del self.executable_dir
        del self.test_data_dir
        del self.sample_name
        del self.peptide_sequence_length

    def test_generate_fasta_runs_and_produces_expected_output(self):
        generate_variant_sequences_input_file  = os.path.join(self.test_data_dir, 'annotated_variants.tsv')
        generate_variant_sequences_output_file = tempfile.NamedTemporaryFile().name
        generate_variant_sequences_executable  = os.path.join(self.executable_dir, 'generate_fasta.py')
        self.assertTrue(py_compile.compile(generate_variant_sequences_executable))

        generate_variant_sequences_command = "%s %s %s %s %s" % (self.python, generate_variant_sequences_executable, generate_variant_sequences_input_file, self.peptide_sequence_length, generate_variant_sequences_output_file)

        self.assertFalse(call(generate_variant_sequences_command, shell=True))
        expected_output_file = os.path.join(self.test_data_dir, ("%s_%s.fa" % (self.sample_name, self.peptide_sequence_length)))
        self.assertTrue(cmp(generate_variant_sequences_output_file, expected_output_file))

if __name__ == '__main__':
    unittest.main()
