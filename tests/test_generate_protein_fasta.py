import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
from tools.pvacseq import generate_protein_fasta

class GenerateFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'tools', 'pvacseq')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_protein_fasta.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'generate_protein_fasta')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_file_with_peptide_sequence_length_17_generates_expected_file(self):
        peptide_sequence_length            = '21'
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            peptide_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))
        expected_tsv_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(generate_protein_fasta_output_tsv, expected_tsv_file))

    def test_output_peptide_sequence_length_longer_that_wildtype(self):
        peptide_sequence_length            = '600'
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            peptide_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
        ], shell=False))
