import unittest
import os
import sys
import tempfile
from filecmp import cmp
import py_compile
import logging
from testfixtures import LogCapture, StringComparison as S

from pvactools.lib.input_file_converter import VcfConverter, FusionInputConverter
from .test_utils import *

class InputFileConverterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'input_file_converter.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'input_file_converter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_truncated_vcf_middle_of_entry(self):
        with self.assertRaises(Exception) as context:
            convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_truncated_middle.vcf')
            convert_vcf_output_file = tempfile.NamedTemporaryFile()

            convert_vcf_params = {
                'input_file'        : convert_vcf_input_file,
                'output_file'       : convert_vcf_output_file.name,
            }
            converter = VcfConverter(**convert_vcf_params)
            converter.execute()
            self.assertTrue("VCF is truncated in the middle of an entry near string " in str(context.exception))

    def test_error_truncated_vcf_end_of_file(self):
        with self.assertRaises(Exception) as context:
            convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_truncated_end.vcf')
            convert_vcf_output_file = tempfile.NamedTemporaryFile()

            convert_vcf_params = {
                'input_file'        : convert_vcf_input_file,
                'output_file'       : convert_vcf_output_file.name,
            }
            converter = VcfConverter(**convert_vcf_params)
            converter.execute()
            self.assertTrue('VCF is truncated at the end of the file' in str(context.exception))

    def test_wildtype_protein_sequence_with_stop(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_wildtype_sequence_with_stop.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()
        logging.disable(logging.NOTSET)

        with LogCapture() as l:
            convert_vcf_params = {
                'input_file'        : convert_vcf_input_file,
                'output_file'       : convert_vcf_output_file.name,
            }
            converter = VcfConverter(**convert_vcf_params)
            converter.execute()
            warn_message = "Transcript WildtypeProtein sequence contains internal stop codon. These can occur in Ensembl transcripts of the biotype polymorphic_pseudogene. Skipping."
            logged_str = "".join(l.actual()[0])
            self.assertTrue(warn_message in logged_str)
            expected_output_file = os.path.join(self.test_data_dir, 'output_wildtype_sequence_with_stop.tsv')
            self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_tx_annotation_generates_expected_tsv(self):
        convert_vcf_input_file              = os.path.join(self.test_data_dir, 'input.tx.vcf')
        convert_vcf_output_file             = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_tx.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_with_empty_tx_and_gx_annotations_generates_expected_tsv(self):
        convert_vcf_input_file              = os.path.join(self.test_data_dir, 'input.empty_tx_gx.vcf')
        convert_vcf_output_file             = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_empty_tx_gx.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_gx_annotation_generates_expected_tsv(self):
        convert_vcf_input_file              = os.path.join(self.test_data_dir, 'input.gx.vcf')
        convert_vcf_output_file             = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_gx.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_multiple_transcripts_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_multiple_transcripts_per_alt_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts_per_alt.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts_per_alt.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_beginning_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_beginning_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_beginning_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_end_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_end_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_end_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_position_out_of_bounds_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_position_out_of_bounds.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_position_out_of_bounds.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_short_wildtype_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_short_wildtype_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_wildtype_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_elongation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_elongation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_truncation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_insertion_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_insertion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_insertion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_deletion_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_deletion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_deletion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_conflicting_alts_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_conflicting_alts.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_conflicting_alts.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_deletion_and_dash_csq_allele(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_dash_csq_allele.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_dash_csq_allele.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_uncalled_genotype_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_uncalled_genotype.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_uncalled_genotype.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_hom_ref_genotype_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_hom_ref_genotype.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_hom_ref_genotype.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_duplicate_index(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_duplicate_index.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_duplicate_index.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_readcount_tags(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input.readcount.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
            'sample_name'       : 'H_NJ-HCC1395-HCC1395',
            'normal_sample_name': 'H_NJ-HCC1395-HCC1396',
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_readcounts.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_missing_csq_format_field_for_variant(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input.no_csq.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_no_csq.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_tsl_vep_field(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_tsl.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_tsl.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_sv_record(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_sv.vcf.gz')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_sv.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_total_length_protein_position(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_total_length.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
            'sample_name'       : 'TUMOR',
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_total_length.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_agfusion_input_generates_expected_tsv(self):
        convert_input_file  = os.path.join(self.test_data_dir, 'agfusion')
        convert_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_input_file,
            'output_file'                : convert_output_file.name,
        }
        converter = FusionInputConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.tsv')
        self.assertTrue(compare(convert_output_file.name, expected_output_file))

    def test_proximal_variants_input(self):
        convert_input_file = os.path.join(self.test_data_dir, 'somatic.vcf.gz')
        convert_input_proximal_variants_file = os.path.join(self.test_data_dir, 'phased.vcf.gz')
        convert_output_file = tempfile.NamedTemporaryFile()
        convert_output_proximal_variants_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file': convert_input_file,
            'output_file': convert_output_file.name,
            'proximal_variants_vcf': convert_input_proximal_variants_file,
            'proximal_variants_tsv': convert_output_proximal_variants_file.name,
            'flanking_bases': 90,
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_proximal_variants_tsv = os.path.join(self.test_data_dir, 'output_proximal_variants.tsv')
        self.assertTrue(cmp(convert_output_proximal_variants_file.name, expected_proximal_variants_tsv))

    def test_protein_altering_variants(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_protein_altering_variants.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'        : convert_vcf_input_file,
            'output_file'       : convert_vcf_output_file.name,
            'sample_name'       : 'TUMOR',
        }
        converter = VcfConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_protein_altering_variants.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))
