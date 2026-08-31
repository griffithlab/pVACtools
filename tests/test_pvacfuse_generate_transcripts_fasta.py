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

from pvactools.lib.generate_transcripts_fasta import PvacfuseGenerateTranscriptsFasta
from pvactools.tools.pvacfuse import generate_transcripts_fasta
from tests.utils import *

class PvacfuseGenerateTranscriptsFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacfuse')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_transcripts_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacfuse_generate_transcripts_fasta')
        cls.transcript_fasta = os.path.join(cls.test_data_dir, 'Homo_sapiens.GRCh38.95.cds.all.fa.gz')

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
        self.assertFalse(result.returncode, "Failed `pvacfuse generate_transcripts_fasta -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_generate_transcripts_fasta_runs(self):
        input_file = os.path.join(self.test_data_dir, 'agfusion')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(generate_transcripts_fasta.main([input_file, self.transcript_fasta, output_file.name]))

    def test_agfusion_input_file_generates_expected_file(self):
        input_file  = os.path.join(self.test_data_dir, 'agfusion')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input': input_file,
            'ref_fasta': self.transcript_fasta,
            'output_file': output_file.name,
            'downstream_sequence_length': None,
        }
        generator = PvacfuseGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_arriba_input_file_generates_expected_file(self):
        input_file  = os.path.join(self.test_data_dir, 'arriba.tsv')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input': input_file,
            'ref_fasta': self.transcript_fasta,
            'output_file': output_file.name,
            'downstream_sequence_length': None,
        }
        generator = PvacfuseGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_arriba.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_downstream_sequence_length(self):
        input_file  = os.path.join(self.test_data_dir, 'agfusion')
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input': input_file,
            'ref_fasta': self.transcript_fasta,
            'output_file': output_file.name,
            'downstream_sequence_length': 50,
        }
        generator = PvacfuseGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.downstream_sequence_length.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))
