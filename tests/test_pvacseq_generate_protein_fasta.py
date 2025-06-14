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

from pvactools.tools.pvacseq import generate_protein_fasta
from tests.utils import *

class GenerateFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'tools', 'pvacseq')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_protein_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'pvacseq_generate_protein_fasta')
        cls.flanking_sequence_length = '10'

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
        self.assertFalse(result.returncode, "Failed `pvacseq generate_protein_fasta -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_generate_protein_fasta_runs(self):
        input_file = os.path.join(self.test_data_dir, 'input.vcf')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(generate_protein_fasta.main([input_file, "25", output_file.name]))
        os.unlink("{}.manufacturability.tsv".format(output_file.name))

    def test_input_vcf_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '--biotypes', 'protein_coding,IG_V_gene',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))
        expected_tsv_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(generate_protein_fasta_output_tsv, expected_tsv_file))
        os.unlink(generate_protein_fasta_output_tsv)

    def test_input_vcf_multi_sample_generates_expected_file(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input_multi_sample.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta_output_tsv = "{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name)

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '-s', 'H_NJ-HCC1395-HCC1395',
            '--biotypes', 'protein_coding,IG_V_gene',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))
        expected_tsv_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(generate_protein_fasta_output_tsv, expected_tsv_file))
        os.unlink(generate_protein_fasta_output_tsv)

    def test_mutant_only(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '--mutant-only',
            '--biotypes', 'protein_coding,IG_V_gene',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutant_only.fasta')
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

    def test_input_tsv(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
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

    def test_input_aggregated_tsv(self):
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_input_tsv   = os.path.join(self.test_data_dir, 'input.aggregated.tsv')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            self.flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '--input-tsv', generate_protein_fasta_input_tsv,
            '--aggregate-report-evaluation', 'Accept',
            '--aggregate-report-evaluation', 'Pending',
            '--biotypes', 'protein_coding,IG_V_gene',
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_with_aggregated_tsv.fasta')
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

    def test_phase_proximal_variants_vcf(self):
        generate_protein_fasta_input_file = os.path.join(self.test_data_dir, 'input_somatic.vcf.gz')
        generate_protein_fasta_phased_proximal_variants_vcf = os.path.join(self.test_data_dir, 'phased.vcf.gz')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            '7',
            generate_protein_fasta_output_file.name,
            '-d', 'full',
            '--phased-proximal-variants-vcf', generate_protein_fasta_phased_proximal_variants_vcf,
        ], shell=False))
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
        proximal_output_file = "{}.proximal_variants.tsv".format(generate_protein_fasta_output_file.name)
        expected_proximal_output_file = os.path.join(self.test_data_dir, 'output_with_phased_vcf.proximal_variants.tsv')
        self.assertTrue(cmp(proximal_output_file, expected_proximal_output_file))
        os.unlink(proximal_output_file)
        expected_output_file = os.path.join(self.test_data_dir, 'output_with_phased_vcf.fasta')
        self.assertTrue(cmp(generate_protein_fasta_output_file.name, expected_output_file))

    def test_output_peptide_sequence_length_longer_that_wildtype(self):
        flanking_sequence_length           = '300'
        generate_protein_fasta_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        generate_protein_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            generate_protein_fasta_input_file,
            flanking_sequence_length,
            generate_protein_fasta_output_file.name,
            '-d', 'full',
        ], shell=False))
        os.unlink("{}.manufacturability.tsv".format(generate_protein_fasta_output_file.name))
