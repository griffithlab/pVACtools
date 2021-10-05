import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

from pvactools.tools.pvacfuse import generate_protein_fasta
from .test_utils import *

class GenerateFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacfuse')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_protein_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacfuse_generate_protein_fasta')
        cls.flanking_sequence_length = '10'

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_agfusion_input_file_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'agfusion')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.fasta')
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

    def test_input_tsv(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'agfusion')
        generate_protein_fasta_input_tsv   = os.path.join(self.test_data_dir, 'input.tsv')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '--input-tsv', generate_protein_fasta_input_tsv,
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_with_tsv.fasta')
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))
