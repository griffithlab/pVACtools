import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
from subprocess import run as subprocess_run
from subprocess import PIPE
import re

from pvactools.lib.generate_protein_fasta import PvacfuseGenerateProteinFasta
from pvactools.tools.pvacfuse import generate_protein_fasta
from tests.utils import *

class PvacfuseGenerateProteinFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacfuse')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_protein_fasta.py')
        cls.test_input_data_dir = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacfuse_generate_transcripts_fasta')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacfuse_generate_protein_fasta')
        cls.flanking_sequence_length = 10

    def test_command(self):
        pvac_script_path = os.path.join(
            self.executable_dir,
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'generate_protein_fasta',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacfuse generate_protein_fasta -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_generate_protein_fasta_runs(self):
        input_file = os.path.join(self.test_input_data_dir, 'output_agfusion.fasta')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(generate_protein_fasta.main([input_file, "25", output_file.name]))
        os.unlink("{}.manufacturability.tsv".format(output_file.name))

    def test_agfusion_input_file_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, 'output_agfusion.fasta')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        params = {
            'transcripts_fasta': generate_protein_fasta_input_file,
            'flanking_sequence_length': self.flanking_sequence_length,
            'output_file': generate_protein_fasta_output_file.name,
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))

    def test_arriba_input_file_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, 'output_arriba.fasta')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        params = {
            'transcripts_fasta': generate_protein_fasta_input_file,
            'flanking_sequence_length': self.flanking_sequence_length,
            'output_file': generate_protein_fasta_output_file.name,
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_arriba.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))

    def test_input_tsv_full(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, 'output_agfusion.fasta')
        generate_protein_fasta_input_tsv   = os.path.join(self.test_data_dir, 'input.tsv')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        params = {
            'transcripts_fasta': generate_protein_fasta_input_file,
            'flanking_sequence_length': self.flanking_sequence_length,
            'output_file': generate_protein_fasta_output_file.name,
            'input_tsv': generate_protein_fasta_input_tsv,
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_with_tsv.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))

    def test_input_tsv_aggregated(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, 'output_agfusion.fasta')
        generate_protein_fasta_input_tsv   = os.path.join(self.test_data_dir, 'input.aggregated.tsv')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        params = {
            'transcripts_fasta': generate_protein_fasta_input_file,
            'flanking_sequence_length': self.flanking_sequence_length,
            'output_file': generate_protein_fasta_output_file.name,
            'input_tsv': generate_protein_fasta_input_tsv,
            'aggregate_report_evaluation': ['Accept', 'Pending'],
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_with_aggregated_tsv.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))

    def test_mutant_only(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, 'output_agfusion.fasta')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        params = {
            'transcripts_fasta': generate_protein_fasta_input_file,
            'flanking_sequence_length': self.flanking_sequence_length,
            'output_file': generate_protein_fasta_output_file.name,
            'mutant_only': True
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.mutant_only.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
