import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

from pvactools.tools.pvacsplice import generate_protein_fasta
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
        'pvacsplice_generate_protein_fasta'
    )

class GenerateFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir       = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacsplice')
        cls.executable           = os.path.join(cls.executable_dir, 'generate_protein_fasta.py')
        cls.test_input_data_dir  = test_input_data_directory()
        cls.test_output_data_dir = test_output_data_directory()
        cls.flanking_sequence_length = '10'

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_vcf_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        generate_protein_fasta_input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        generate_protein_fasta_input_fasta = os.path.join(self.test_input_data_dir, "inputs", "all_sequences_chr1.fa.gz")
        unzipped_fasta_file = gunzip_file(generate_protein_fasta_input_fasta)
        generate_protein_fasta_input_gtf   = os.path.join(self.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz")
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv  = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            generate_protein_fasta_input_vcf,
            unzipped_fasta_file,
            generate_protein_fasta_input_gtf,
            '-s', 'HCC1395_TUMOR_DNA',
        ], shell=False))
        expected_output_file = os.path.join(self.test_output_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        expected_tsv_file = os.path.join(self.test_output_data_dir, 'output.tsv')
        self.assertTrue(cmp(generate_protein_fasta_output_tsv, expected_tsv_file))
        os.unlink(generate_protein_fasta_output_tsv)
        os.unlink(unzipped_fasta_file)

    def test_input_tsv(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        generate_protein_fasta_input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        generate_protein_fasta_input_fasta = os.path.join(self.test_input_data_dir, "inputs", "all_sequences_chr1.fa.gz")
        unzipped_fasta_file = gunzip_file(generate_protein_fasta_input_fasta)
        generate_protein_fasta_input_gtf   = os.path.join(self.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz")
        generate_protein_fasta_input_tsv   = os.path.join(self.test_output_data_dir, "Test.filtered.tsv")
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv  = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            generate_protein_fasta_input_vcf,
            unzipped_fasta_file,
            generate_protein_fasta_input_gtf,
            '-s', 'HCC1395_TUMOR_DNA',
            '--input-tsv', generate_protein_fasta_input_tsv,
        ], shell=False))
        expected_output_file = os.path.join(self.test_output_data_dir, 'output_with_tsv.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        os.unlink(unzipped_fasta_file)

    def test_input_aggregated_tsv(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_input_data_dir, "inputs", "splice_junctions_chr1.tsv")
        generate_protein_fasta_input_vcf   = os.path.join(self.test_input_data_dir, "inputs", "annotated.expression_chr1.vcf.gz")
        generate_protein_fasta_input_fasta = os.path.join(self.test_input_data_dir, "inputs", "all_sequences_chr1.fa.gz")
        unzipped_fasta_file = gunzip_file(generate_protein_fasta_input_fasta)
        generate_protein_fasta_input_gtf   = os.path.join(self.test_input_data_dir, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz")
        generate_protein_fasta_input_tsv   = os.path.join(self.test_output_data_dir, "Test.all_epitopes.aggregated.tsv")
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv  = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            generate_protein_fasta_input_vcf,
            unzipped_fasta_file,
            generate_protein_fasta_input_gtf,
            '-s', 'HCC1395_TUMOR_DNA',
            '--input-tsv', generate_protein_fasta_input_tsv,
            '--aggregate-report-evaluation', 'Pending'
        ], shell=False))
        expected_output_file = os.path.join(self.test_output_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        os.unlink(unzipped_fasta_file)
