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

from pvactools.tools.pvacsplice import generate_transcripts_fasta
from pvactools.lib.generate_transcripts_fasta import PvacspliceGenerateTranscriptsFasta
from tests.utils import *

def test_input_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacsplice'
    )

def test_output_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacsplice_generate_transcripts_fasta'
    )

class PvacspliceGenerateTranscriptsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir       = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacsplice')
        cls.executable           = os.path.join(cls.executable_dir, 'generate_transcripts_fasta.py')
        cls.test_input_data_dir  = test_input_data_directory()
        cls.test_output_data_dir = test_output_data_directory()
        cls.gtf_file_chr1 = os.path.join(cls.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.filtered.gtf")
        cls.gtf_file_chr2 = os.path.join(cls.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105.chr2.filtered.gtf")
        cls.gtf_file_chr3 = os.path.join(cls.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105.chr3.filtered.gtf")
        cls.fasta_file_chr1 = gunzip_file(os.path.join(cls.test_input_data_dir, "inputs", "all_sequences_chr1.fa.gz"))
        cls.fasta_file_chr2 = gunzip_file(os.path.join(cls.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz"))
        cls.fasta_file_chr3 = gunzip_file(os.path.join(cls.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz"))

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.fasta_file_chr1)
        os.unlink(cls.fasta_file_chr2)
        os.unlink(cls.fasta_file_chr3)

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
        self.assertFalse(result.returncode, "Failed `pvacsplice generate_transcripts_fasta -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_generate_transcripts_fasta_runs(self):
        input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(generate_transcripts_fasta.main([
            input_file,
            output_file.name,
            input_vcf,
            self.fasta_file_chr1,
            self.gtf_file_chr1,
            '-s', 'HCC1395_TUMOR_DNA',
        ]))

    def test_input_vcf_generates_expected_file(self):
        input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_file': input_file,
            'output_file': output_file.name,
            'annotated_vcf': input_vcf,
            'ref_fasta': self.fasta_file_chr1,
            'gtf_file': self.gtf_file_chr1,
            'sample_name': 'HCC1395_TUMOR_DNA',
        }
        generator = PvacspliceGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_output_data_dir, 'output.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_input_short_sequence_generates_expected_file(self):
        input_file  = os.path.join(self.test_input_data_dir, "inputs", "regtools.short_sequence.tsv")
        input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "short_sequence.vcf.gz")
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_file': input_file,
            'output_file': output_file.name,
            'annotated_vcf': input_vcf,
            'ref_fasta': self.fasta_file_chr3,
            'gtf_file': self.gtf_file_chr3,
            'sample_name': 'TumorDNA',
        }
        generator = PvacspliceGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_output_data_dir, 'output.short.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_input_unsupported_amino_acid_generates_expected_file(self):
        input_file  = os.path.join(self.test_input_data_dir, "inputs", "regtools.unsupported_aa.tsv")
        input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "unsupported_aa.vcf.gz")
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_file': input_file,
            'output_file': output_file.name,
            'annotated_vcf': input_vcf,
            'ref_fasta': self.fasta_file_chr2,
            'gtf_file': self.gtf_file_chr2,
            'sample_name': 'TumorDNA',
        }
        generator = PvacspliceGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_output_data_dir, 'output.unsupported_aa.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))

    def test_downstream_sequence_length_generates_expected_file(self):
        input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        output_file = tempfile.NamedTemporaryFile()

        params = {
            'input_file': input_file,
            'output_file': output_file.name,
            'annotated_vcf': input_vcf,
            'ref_fasta': self.fasta_file_chr1,
            'gtf_file': self.gtf_file_chr1,
            'sample_name': 'HCC1395_TUMOR_DNA',
            'downstream_sequence_length': 10,
        }
        generator = PvacspliceGenerateTranscriptsFasta(**params)
        generator.execute()

        expected_output_file = os.path.join(self.test_output_data_dir, 'output.downstream_sequence_length.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))
