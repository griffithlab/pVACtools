import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class ConvertVcfTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'convert_vcf.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'convert_vcf')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_vcf_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutliple_transcripts_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_multiple_transcripts_per_alt_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts_per_alt.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts_per_alt.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_beginning_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_beginning_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_beginning_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_end_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_end_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_end_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_position_out_of_bounds_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_position_out_of_bounds.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_position_out_of_bounds.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_short_wildtype_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_short_wildtype_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_wildtype_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_elongation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_elongation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_truncation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
            ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_insertion_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_insertion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_insertion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_deletion_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_deletion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, convert_vcf_input_file, convert_vcf_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_deletion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))
