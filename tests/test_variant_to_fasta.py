import sys
import os
import unittest
import tempfile
from filecmp import cmp
import py_compile
import itertools

from pvactools.lib.variant_to_fasta import VariantToFasta
from tests.utils import *

class VariantToFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'variant_to_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'variant_to_fasta')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_file_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_mutation_at_relative_end_of_full_sequence_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_mutation_at_relative_end_of_full_sequence.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_end_of_full_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_mutation_at_relative_beginning_of_full_sequence_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_mutation_at_relative_beginning_of_full_sequence.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_beginning_of_full_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_wildtype_sequence_shorter_than_desired_peptide_sequence_length_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_short_wildtype_sequence.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_wildtype_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_position_out_of_bounds_generates_empty_output_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_position_out_of_bounds.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.stat(generate_fasta_output_file.name).st_size, 0)

    def test_input_file_with_inframe_insertion_amino_acid_replacement_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_replacement.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_deletion_amino_acid_replacement_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_replacement.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_replacement.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_insertion_amino_acid_insertion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_insertion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_deletion_amino_acid_deletion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_deletion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_deletion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_inframe_deletion_range(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_deletion_range.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_range.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_feature_truncation_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_feature_truncation2_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation2.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': 100,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation2.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_feature_elongation_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_elongation.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_range_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_range.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_range.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_resulting_short_fasta_sequence(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_short_fasta_sequence.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_fasta_sequence.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_downstream_sequence_length_limit_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': 20,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_downstream_sequence_length_limit.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_dnp_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_dnp.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_dnp.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_protein_change_with_asterisk_in_wildtype(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'protein_change_with_asterisk_in_wildtype')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.path.getsize(generate_fasta_output_file.name), 0)

    def test_proximal_variants_generate_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_somatic_variant_with_proximal_variants.tsv')
        generate_fasta_proximal_variants_file = os.path.join(self.test_data_dir, 'input_proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_proximal_variants.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_protein_change_with_asterisk_in_wildtype_and_mutant(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'protein_change_with_asterisk_in_wildtype_and_mutant')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.path.getsize(generate_fasta_output_file.name), 0)

    def test_protein_change_with_X_in_wildtype_and_mutatnt(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'protein_change_with_X_in_wildtype_and_mutant')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.path.getsize(generate_fasta_output_file.name), 0)

    def test_proximal_variants_for_inframe_insertion_generate_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.tsv')
        generate_fasta_proximal_variants_file = os.path.join(self.test_data_dir, 'input_proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_proximal_variants_inframe_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_protein_change_with_multiple_asterisks(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'protein_change_with_multiple_asterisks')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.path.getsize(generate_fasta_output_file.name), 0)

    def test_inframe_insertion_with_no_aa_change(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_no_aa_change.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        self.assertEqual(os.path.getsize(generate_fasta_output_file.name), 0)

    def test_proximal_variant_on_same_codon_as_somatic_variant_results_in_novel_peptide(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'proximal_variant_on_same_codon_as_somatic_variant_results_in_novel_peptide')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_proximal_variants_file = os.path.join(test_data_dir, 'proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_proximal_variant_on_same_codon_as_somatic_variant_results_in_stop_codon(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'proximal_variant_on_same_codon_as_somatic_variant_results_in_stop_codon')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_proximal_variants_file = os.path.join(test_data_dir, 'proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(test_data_dir, 'output_proximal_variant_on_same_codon_as_somatic_variant_results_in_stop_codon.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_multiple_proximal_variants_on_same_codon_results_in_novel_peptide(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'multiple_proximal_variants_on_same_codon_results_in_novel_peptide')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_proximal_variants_file = os.path.join(test_data_dir, 'proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_multiple_proximal_variants_on_same_codon_results_in_stop_codon(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'multiple_proximal_variants_on_same_codon_results_in_stop_codon')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_proximal_variants_file = os.path.join(test_data_dir, 'proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(test_data_dir, 'output_multiple_proximal_variants_on_same_codon_results_in_stop_codon.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_proximal_dnp(self):
        test_data_dir                  = os.path.join(self.test_data_dir, 'proximal_dnp')
        generate_fasta_input_file      = os.path.join(test_data_dir, 'input.tsv')
        generate_fasta_proximal_variants_file = os.path.join(test_data_dir, 'proximal_variants.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'proximal_variants_file'    : generate_fasta_proximal_variants_file,
            'downstream_sequence_length': None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(test_data_dir, 'output.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_protein_altering_variant_insertion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_protein_altering_variant_insertion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_protein_altering_variant_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_very_long_insertion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_very_long_insertion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_very_long_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_complex_inframe_insertion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_complex_inframe_insertion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_complex_inframe_insertion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_frameshift_variant_position_1_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_frameshift_variant_position_1.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_position_1.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

    def test_input_file_with_repetitive_deletion_generates_expected_file(self):
        generate_fasta_input_file      = os.path.join(self.test_data_dir, 'input_repetitive_deletion.tsv')
        generate_fasta_output_file     = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'                : generate_fasta_input_file,
            'output_file'               : generate_fasta_output_file.name,
            'downstream_sequence_length': None,
            'proximal_variants_file'    : None,
        }
        generator = VariantToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_repetitive_deletion.fasta')
        self.assertTrue(cmp(generate_fasta_output_file.name, expected_output_file))

if __name__ == '__main__':
    unittest.main()
