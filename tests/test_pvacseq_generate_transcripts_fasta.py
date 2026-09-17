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
import shutil

from pvactools.tools.pvacseq import generate_transcripts_fasta
from pvactools.lib.generate_transcripts_fasta import PvacseqGenerateTranscriptsFasta
from tests.utils import *

class PvacseqGenerateTranscriptsFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacseq')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_transcripts_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacseq_generate_transcripts_fasta')
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
            'generate_transcripts_fasta',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq generate_transcripts_fasta -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_generate_transcripts_fasta_runs(self):
        input_file = os.path.join(self.test_data_dir, 'input.vcf')
        output_file = tempfile.NamedTemporaryFile()
        generate_transcripts_fasta.main([input_file, output_file.name])

    def test_input_vcf_generates_expected_file(self):
        input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_vcf': input_file,
            'output_file': output_file.name,
            'downstream_sequence_length': None,
            'biotypes': ['protein_coding', 'IG_V_gene'],
        }
        generator = PvacseqGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_input_vcf_multi_sample_generates_expected_file(self):
        input_file  = os.path.join(self.test_data_dir, 'input_multi_sample.vcf')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_vcf': input_file,
            'output_file': output_file.name,
            'downstream_sequence_length': None,
            'sample_name': 'H_NJ-HCC1395-HCC1395',
            'biotypes': ['protein_coding', 'IG_V_gene'],
        }
        generator = PvacseqGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_phase_proximal_variants_vcf(self):
        input_file = os.path.join(self.test_data_dir, 'input_somatic.vcf.gz')
        phased_proximal_variants_vcf = os.path.join(self.test_data_dir, 'phased.vcf.gz')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_vcf': input_file,
            'phased_proximal_variants_vcf': phased_proximal_variants_vcf,
            'allow_incomplete_transcripts': True,
            'output_file': output_file.name,
        }
        generator = PvacseqGenerateTranscriptsFasta(**params)
        generator.generate_fasta()
        shutil.copy(generator.fasta_file_path, output_file.name)

        expected_output_file = os.path.join(self.test_data_dir, 'output_with_phased_vcf.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))
        expected_output_file = os.path.join(self.test_data_dir, 'output_with_phased_vcf.proximal_variants.tsv')
        self.assertTrue(cmp(os.path.join(generator.temp_dir, 'tmp.proximal_variants.tsv'), expected_output_file))
