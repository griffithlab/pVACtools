import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
from pvacseq.lib import generate_fasta

class GenerateFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'generate_fasta.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'generate_fasta')
        cls.peptide_sequence_length = '21'

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_file_with_peptide_sequence_length_17_generates_expected_file(self):
        peptide_sequence_length = '17'
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_peptide_sequence_length_17.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_peptide_sequence_length_21_generates_expected_file(self):
        peptide_sequence_length = '21'
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_peptide_sequence_length_21.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_peptide_sequence_length_31_generates_expected_file(self):
        peptide_sequence_length = '31'
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_peptide_sequence_length_31.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_mutation_at_relative_end_of_full_sequence_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_end_of_full_sequence.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_end_of_full_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_mutation_at_relative_beginning_of_full_sequence_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_beginning_of_full_sequence.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_beginning_of_full_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_wildtype_sequence_shorter_than_desired_peptide_sequence_length_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_short_wildtype_sequence.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_wildtype_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_multiple_transcripts_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_multiple_transcripts_per_alternate_allele_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts_per_alt.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts_per_alt.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_position_out_of_bounds_generates_empty_output_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_position_out_of_bounds.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        self.assertEqual(os.stat(generate_fasta_output_file.name).st_size, 0)

    def test_input_file_with_inframe_insertion_amino_acid_replacement_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_replacement.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_deletion_amino_acid_replacement_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_replacement.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_replacement.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_insertion_amino_acid_insertion_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_insertion.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_deletion_amino_acid_deletion_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_deletion.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_deletion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_feature_truncation_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_feature_elongation_generates_expected_file(self):
        generate_fasta_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.tsv')
        generate_fasta_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, generate_fasta_input_file, self.peptide_sequence_length, generate_fasta_output_file.name
        ], shell=False))
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_elongation.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_distance_from_start_works_as_expected(self):
        sequence = 'KKLKILGMPFRNIRSILKMVN'
        position = 5
        self.assertEqual(generate_fasta.distance_from_start(position, sequence), 5)

    def test_distance_from_end_works_as_expected(self):
        sequence = 'KKLKILGMPFRNIRSILKMVN'
        position = 5
        self.assertEqual(generate_fasta.distance_from_end(position, sequence), 15)

if __name__ == '__main__':
    unittest.main()
